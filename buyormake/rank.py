"""Cost-aware ranking and proposal filtering for candidate molecules."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import httpx

from .cache import DEFAULT_CACHE_PATH, DEFAULT_ONLINE_TTL_SECONDS, ResponseCache
from .canonical import normalize_smiles
from .cost import (
    cheapest_offer_for_target,
    estimate_route_cost_bounds,
    min_outlay_for_offers,
    skipped_leaves_from_method,
)
from .handoff import _min_offer_lead_days
from .routes import (
    DEFAULT_REACTION_DB,
    AiZynthFinderBackend,
    ReactionLookupBackend,
    TemplateFallbackBackend,
    find_routes,
)
from .schema import (
    AvailabilityRecord,
    FilterSummary,
    InChIKeyMatchLiteral,
    Molecule,
    Offer,
    PriceStatusLiteral,
    ProposalDecision,
    Provenance,
    Quantity,
    RankResult,
    Route,
    SAScores,
    SoftwareRef,
)
from .scores import DEFAULT_SCSCORE_WEIGHTS, compute_scores
from .sources import SOURCE_BY_ID
from .stock import LocalStockBackend, PubChemVendorBackend

_DEFAULT_STOCK_FILES = (
    Path("data/stock.sqlite"),
    Path("data/stock.tsv"),
    Path("data/stock.jsonl"),
)

_STEP_SOURCE_TO_SOURCE_ID = {
    "ord_lookup": "ord_data",
    "uspto_lookup": "uspto_lowe",
    "aizynthfinder": "aizynthfinder_models",
    "paroutes_lookup": "paroutes",
    "orderly_lookup": "orderly",
    "syntheseus_uspto_50k_lookup": "syntheseus_uspto_50k",
    "synplanner_lookup": "synplanner_data",
}


def _rank_sort_key(res: RankResult) -> tuple[int, int, int, str, float]:
    if res.status == "buyable":
        if res.best_price is not None and res.best_price.value is not None:
            price_val = float(res.best_price.value)
            unit = res.best_price.unit
            if unit == "USD/g":
                return (0, 0, 0, unit, price_val)
            return (0, 1, 0, unit, price_val)
        return (0, 2, 0, "", 0.0)

    if res.status == "makeable":
        steps = res.n_steps if res.n_steps is not None else 999
        has_cost = res.est_cost is not None and res.est_cost.value is not None
        if has_cost and res.est_cost is not None and res.est_cost.value is not None:
            unit = res.est_cost.unit
            cost_val = float(res.est_cost.value)
            if res.cost_complete and unit in {"USD", "USD/g"}:
                cost_tier = 0
            elif res.cost_complete:
                cost_tier = 1
            else:
                cost_tier = 2
            return (1, steps, cost_tier, unit, cost_val)
        return (1, steps, 3, "", 0.0)

    if res.status == "unknown":
        return (2, 0, 0, "", 0.0)
    return (3, 0, 0, "", 0.0)


def _unique_urls(offers: Iterable[Offer]) -> list[str]:
    seen: set[str] = set()
    urls: list[str] = []
    for off in offers:
        if off.source_url and off.source_url not in seen:
            seen.add(off.source_url)
            urls.append(off.source_url)
    return urls


def _collect_provenance_and_licenses(
    offers: list[Offer],
    routes: list[Route],
    sa_scores: SAScores,
    offer_cache: dict[str, list[Offer]],
) -> tuple[list[str], list[str]]:
    """Aggregate all contributing data/model `source_id`s and licenses (target, route steps, and leaves)."""
    used_set: set[str] = {o.source_id for o in offers if o.source_id}
    lic_set: set[str] = {o.license for o in offers if o.license}

    for r in routes:
        for st in r.steps:
            mapped_src = _STEP_SOURCE_TO_SOURCE_ID.get(st.source) or (
                st.source if st.source in SOURCE_BY_ID else None
            )
            if mapped_src:
                used_set.add(mapped_src)
        for leaf_smi in r.leaves:
            lnm = normalize_smiles(leaf_smi)
            if lnm is None:
                continue
            for lo in offer_cache.get(lnm.inchikey, []):
                if lo.source_id:
                    used_set.add(lo.source_id)
                if lo.license:
                    lic_set.add(lo.license)

    for sid in used_set:
        if sid in SOURCE_BY_ID:
            lic_set.add(SOURCE_BY_ID[sid].license)

    if sa_scores.sascore is not None:
        lic_set.add("BSD-3-Clause")
    if sa_scores.scscore is not None:
        lic_set.add("MIT")
    if sa_scores.rascore is not None:
        lic_set.add("MIT")

    return sorted(used_set), sorted(lic_set)


@dataclass
class _EvaluatedEntry:
    canonical_smiles: str
    inchikey: str
    formula: str
    status: str
    n_steps: Optional[int]
    best_price: Optional[Quantity]
    est_cost: Optional[Quantity]
    est_cost_yield_adjusted: Optional[Quantity]
    min_outlay_cost: Optional[Quantity]
    step_yield_assumed: float
    cost_complete: bool
    lead_time_days: Optional[int]
    moq_pack_g: Optional[float]
    precursor_smiles: list[str]
    inchikey_match: InChIKeyMatchLiteral
    matched_offer_inchikeys: list[str]
    stereoisomer_inchikeys_in_stock: list[str]
    stereoisomer_offers: list[Offer]
    leaf_offers: dict[str, list[Offer]]
    sa_scores: SAScores
    offers: list[Offer]
    routes: list[Route]
    reasons: list[str]
    evidence: list[str]
    used_sources: list[str]
    licenses: list[str]


def evaluate_candidates(
    candidates: list[str],
    *,
    candidate_names: Optional[list[Optional[str]]] = None,
    max_steps: int = 3,
    backends: tuple[str, ...] = ("local_stock", "pubchem"),
    route_backend: str = "auto",
    offline: bool = False,
    cache_path: Path | str = DEFAULT_CACHE_PATH,
    cache_max_age_seconds: Optional[float] = None,
    refresh_cache: bool = False,
    stock_paths: Optional[list[Path | str]] = None,
    reaction_db_path: Path | str = DEFAULT_REACTION_DB,
    scscore_weights: Path | str = DEFAULT_SCSCORE_WEIGHTS,
    aizynth_config: Optional[Path | str] = None,
    include_scores: bool = True,
    desalt_fallback: bool = True,
    stereo_relaxed: bool = False,
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
    target_scale_g: float = 1.0,
    max_pack_multiple: Optional[float] = None,
    currency: str = "USD",
    http_client: Optional[httpx.Client] = None,
) -> tuple[list[RankResult], list[AvailabilityRecord]]:
    """Evaluate candidates and return `(sorted RankResults, AvailabilityRecords in input order)`."""
    effective_ttl = (
        None
        if offline and cache_max_age_seconds is None
        else (
            cache_max_age_seconds
            if cache_max_age_seconds is not None
            else (None if offline else DEFAULT_ONLINE_TTL_SECONDS)
        )
    )
    cache = ResponseCache(
        cache_path,
        max_age_seconds=effective_ttl,
        refresh=refresh_cache,
        readonly=(offline and not refresh_cache),
    )
    resolved_stock = (
        [Path(p) for p in stock_paths]
        if stock_paths is not None
        else [p for p in _DEFAULT_STOCK_FILES if p.exists()]
    )
    local_stock = LocalStockBackend(resolved_stock) if "local_stock" in backends else None
    pubchem = (
        PubChemVendorBackend(cache, offline=offline, client=http_client)
        if "pubchem" in backends
        else None
    )
    lookup_backend = ReactionLookupBackend(reaction_db_path)
    template_backend = TemplateFallbackBackend()
    aizynth_backend = AiZynthFinderBackend(aizynth_config)

    exact_offer_cache: dict[str, list[Offer]] = {}
    leaf_offer_cache: dict[str, list[Offer]] = {}
    offer_mol_weights: dict[str, float] = (
        local_stock.known_mol_weights() if local_stock is not None else {}
    )

    def lookup_exact_offers(inchikey: str, canonical_smiles: str) -> list[Offer]:
        if inchikey in exact_offer_cache:
            return exact_offer_cache[inchikey]
        found: list[Offer] = []
        if local_stock is not None:
            found.extend(local_stock.lookup(inchikey))
            offer_mol_weights.update(local_stock.known_mol_weights())
        if pubchem is not None:
            found.extend(pubchem.lookup(inchikey, canonical_smiles))
        exact_offer_cache[inchikey] = found
        return found

    def lookup_leaf_offers(inchikey: str, canonical_smiles: str) -> list[Offer]:
        if inchikey in leaf_offer_cache:
            return leaf_offer_cache[inchikey]
        lnm = normalize_smiles(canonical_smiles, desalt=False)
        if lnm is not None:
            offer_mol_weights[lnm.inchikey] = lnm.mol_weight
        found = list(lookup_exact_offers(inchikey, canonical_smiles))
        if not found and desalt_fallback:
            # Neutral retrosynthetic leaf -> pharmaceutical salt/hydrate stock equivalent
            if local_stock is not None and (lnm is None or not lnm.has_permanent_charge):
                found.extend(
                    local_stock.lookup_parent_equivalents(inchikey, exclude_inchikey=inchikey)
                )
                offer_mol_weights.update(local_stock.known_mol_weights())
            # Salt/hydrate leaf -> desalted neutral parent or sibling salt stock
            if not found and lnm is not None and (
                "." in lnm.canonical_smiles
                or "+" in lnm.canonical_smiles
                or "-" in lnm.canonical_smiles
            ):
                p_nm = normalize_smiles(canonical_smiles, desalt=True)
                if p_nm is not None and p_nm.inchikey != inchikey:
                    offer_mol_weights[p_nm.inchikey] = p_nm.mol_weight
                    found = list(lookup_exact_offers(p_nm.inchikey, p_nm.canonical_smiles))
                    if not found and local_stock is not None:
                        found = local_stock.lookup_parent_equivalents(
                            p_nm.inchikey, exclude_inchikey=inchikey
                        )
                        offer_mol_weights.update(local_stock.known_mol_weights())
        leaf_offer_cache[inchikey] = found
        return found

    results: list[RankResult] = []
    records: list[AvailabilityRecord] = []
    memo_by_ik: dict[str, _EvaluatedEntry] = {}

    try:
        for idx, raw in enumerate(candidates):
            cand_name = (
                candidate_names[idx]
                if (candidate_names is not None and idx < len(candidate_names))
                else None
            )
            nm = normalize_smiles(raw, desalt=False)
            if nm is None:
                results.append(
                    RankResult(
                        candidate=raw,
                        name=cand_name,
                        canonical_smiles=None,
                        inchikey=None,
                        status="invalid",
                        n_steps=None,
                        best_price=None,
                        est_cost=None,
                        step_yield_assumed=default_step_yield,
                        cost_complete=False,
                        sa_scores=None,
                        reasons=["invalid SMILES; RDKit could not parse structure"],
                        evidence=[],
                    )
                )
                continue

            offer_mol_weights[nm.inchikey] = nm.mol_weight

            entry = memo_by_ik.get(nm.inchikey)
            if entry is None:
                offers = list(lookup_exact_offers(nm.inchikey, nm.canonical_smiles))
                inchikey_match: InChIKeyMatchLiteral = "exact" if offers else "none"
                desalt_reason: Optional[str] = None
                stereo_reason: Optional[str] = None
                parent_nm = None

                # Symmetric salt-form equivalence (salt -> freebase AND freebase -> salt)
                if not offers and desalt_fallback:
                    if (
                        "." in nm.canonical_smiles
                        or "+" in nm.canonical_smiles
                        or "-" in nm.canonical_smiles
                    ):
                        cand_parent = normalize_smiles(raw, desalt=True)
                        if cand_parent is not None and cand_parent.inchikey != nm.inchikey:
                            parent_nm = cand_parent
                            offer_mol_weights[parent_nm.inchikey] = parent_nm.mol_weight
                            parent_offers = list(
                                lookup_exact_offers(
                                    parent_nm.inchikey, parent_nm.canonical_smiles
                                )
                            )
                            if not parent_offers and local_stock is not None:
                                parent_offers = local_stock.lookup_parent_equivalents(
                                    parent_nm.inchikey, exclude_inchikey=nm.inchikey
                                )
                                offer_mol_weights.update(local_stock.known_mol_weights())
                            if parent_offers:
                                offers = parent_offers
                                inchikey_match = "salt_equivalent"
                                desalt_reason = (
                                    f"matched desalted parent structure {parent_nm.canonical_smiles} "
                                    f"({parent_nm.inchikey}) after exact salt/mixture form had 0 offers"
                                )
                    elif local_stock is not None and not nm.has_permanent_charge:
                        salt_offers = local_stock.lookup_parent_equivalents(
                            nm.inchikey, exclude_inchikey=nm.inchikey
                        )
                        offer_mol_weights.update(local_stock.known_mol_weights())
                        if salt_offers:
                            offers = salt_offers
                            inchikey_match = "salt_equivalent"
                            salt_iks = ", ".join(sorted({o.inchikey for o in salt_offers}))
                            desalt_reason = (
                                f"matched pharmaceutical salt/hydrate stock equivalent ({salt_iks}) "
                                f"for neutral parent {nm.canonical_smiles} ({nm.inchikey})"
                            )

                # 14-char skeletal connectivity block lookup (stereoisomers / racemates)
                excluded_stereo_iks = {
                    nm.inchikey,
                    *([parent_nm.inchikey] if parent_nm is not None else []),
                    *(o.inchikey for o in offers),
                }
                stereo_offers = (
                    local_stock.lookup_connectivity_block(
                        nm.connectivity_block, exclude_inchikeys=excluded_stereo_iks
                    )
                    if local_stock is not None
                    else []
                )
                stereo_iks = sorted({o.inchikey for o in stereo_offers})
                if not offers and stereo_offers:
                    if stereo_relaxed:
                        offers = list(stereo_offers)
                        inchikey_match = "stereoisomer_relaxed"
                        stereo_reason = (
                            f"stereoisomer_relaxed: matched 14-char connectivity block "
                            f"{nm.connectivity_block} against stereoisomer/racemate stock "
                            f"({', '.join(stereo_iks)}); stereochemistry differs from query"
                        )
                    else:
                        stereo_reason = (
                            f"stereoisomer_advisory: exact 27-char InChIKey had 0 offers, "
                            f"but connectivity block {nm.connectivity_block} matched "
                            f"stereoisomer/racemate stock ({', '.join(stereo_iks)}); "
                            f"pass stereo_relaxed=True to accept stereoisomeric stock"
                        )

                sa_scores = (
                    compute_scores(
                        nm.canonical_smiles,
                        scscore_weights=scscore_weights,
                        mol=nm.mol,
                    )
                    if include_scores
                    else SAScores()
                )
                routes: list[Route] = []
                rec_leaf_offers: dict[str, list[Offer]] = {}
                reasons: list[str] = []
                if desalt_reason is not None:
                    reasons.append(desalt_reason)
                if stereo_reason is not None:
                    reasons.append(stereo_reason)
                evidence = _unique_urls(offers)

                est_cost_yield_adjusted: Optional[Quantity] = None
                min_outlay_cost: Optional[Quantity] = None
                lead_time_days: Optional[int] = None
                moq_pack_g: Optional[float] = None
                precursor_smiles: list[str] = []
                matched_offer_iks = sorted({o.inchikey for o in offers})

                if offers:
                    status = "buyable"
                    max_pack_g = (
                        target_scale_g * max_pack_multiple
                        if max_pack_multiple is not None
                        else None
                    )
                    best_price = cheapest_offer_for_target(
                        offers,
                        currency=currency,
                        target_mol_weight=nm.mol_weight,
                        target_inchikey=nm.inchikey,
                        offer_mol_weights=offer_mol_weights,
                        max_pack_grams=max_pack_g,
                    )
                    min_outlay_cost, _, moq_pack_g = min_outlay_for_offers(
                        offers,
                        currency=currency,
                        required_grams=target_scale_g,
                        target_mol_weight=nm.mol_weight,
                        target_inchikey=nm.inchikey,
                        offer_mol_weights=offer_mol_weights,
                    )
                    lead_time_days = _min_offer_lead_days(offers)
                    est_cost = None
                    cost_complete = False
                    n_steps = 0
                    vendors_preview = ", ".join(dict.fromkeys(o.vendor for o in offers[:3]))
                    suffix = f" (+{len(offers) - 3} more offers)" if len(offers) > 3 else ""
                    reasons.append(
                        f"listed by {len(offers)} catalog offer(s): {vendors_preview}{suffix}"
                    )
                    if best_price is None:
                        reasons.append(
                            "listed price unknown (source does not publish open catalog prices)"
                        )
                    else:
                        reasons.append(
                            f"cheapest listed price: {best_price.raw_unit_str or best_price.value}"
                        )
                else:
                    route_target = parent_nm if parent_nm is not None else nm
                    routes = find_routes(
                        route_target.canonical_smiles,
                        route_target.inchikey,
                        max_steps=max_steps,
                        route_backend=route_backend,
                        offer_lookup=lookup_leaf_offers,
                        lookup_backend=lookup_backend,
                        template_backend=template_backend,
                        aizynth_backend=aizynth_backend,
                        offer_cache=leaf_offer_cache,
                        target_mol=route_target.mol,
                        currency=currency,
                        desalt_fallback=desalt_fallback,
                        default_step_yield=default_step_yield,
                        stoich_excess=stoich_excess,
                        target_scale_g=target_scale_g,
                        max_pack_multiple=max_pack_multiple,
                        offer_mol_weights=offer_mol_weights,
                    )
                    for r_item in routes:
                        for leaf_smi in r_item.leaves:
                            lnm = normalize_smiles(leaf_smi)
                            if lnm is not None:
                                l_offs_cached = leaf_offer_cache.get(lnm.inchikey)
                                if l_offs_cached:
                                    rec_leaf_offers[lnm.inchikey] = list(l_offs_cached)

                    solved = [r for r in routes if r.leaves and all(r.leaves_purchasable)]
                    if solved:
                        best_route = solved[0]
                        status = "makeable"
                        n_steps = best_route.n_steps
                        best_price = None
                        est_cost = best_route.est_cost
                        est_cost_yield_adjusted = best_route.est_cost_yield_adjusted
                        min_outlay_cost = best_route.min_outlay_cost
                        cost_complete = best_route.est_cost_complete
                        leaf_urls: list[str] = []
                        leaf_leads: list[int] = []
                        missing_leads = False
                        for leaf_smi in dict.fromkeys(best_route.leaves):
                            lnm = normalize_smiles(leaf_smi)
                            if lnm is not None:
                                precursor_smiles.append(lnm.canonical_smiles)
                                l_offs = leaf_offer_cache.get(lnm.inchikey, [])
                                leaf_urls.extend(_unique_urls(l_offs))
                                l_day = _min_offer_lead_days(l_offs)
                                if l_day is None:
                                    missing_leads = True
                                else:
                                    leaf_leads.append(l_day)
                            else:
                                precursor_smiles.append(leaf_smi)
                                missing_leads = True
                        if not missing_leads and leaf_leads:
                            lead_time_days = max(leaf_leads)
                        breakdown_moq = estimate_route_cost_bounds(
                            route_target.canonical_smiles,
                            best_route.leaves,
                            rec_leaf_offers,
                            currency=currency,
                            steps=best_route.steps,
                            default_step_yield=default_step_yield,
                            stoich_excess=stoich_excess,
                            target_scale_g=target_scale_g,
                            max_pack_multiple=max_pack_multiple,
                            offer_mol_weights=offer_mol_weights,
                        )
                        moq_pack_g = breakdown_moq.moq_pack_g
                        evidence = list(dict.fromkeys(leaf_urls))
                        if parent_nm is not None:
                            reasons.append(
                                f"evaluated retrosynthesis on desalted parent {parent_nm.canonical_smiles} "
                                f"({parent_nm.inchikey})"
                            )
                        reasons.append(
                            f"makeable in {best_route.n_steps} step(s) via {best_route.method} "
                            f"from {len(best_route.leaves)} purchasable leaf/leaves"
                        )
                        if est_cost is None:
                            reasons.append("leaf prices unknown; route cost could not be estimated")
                        elif not cost_complete:
                            reasons.append(
                                f"partial route cost lower bound {est_cost.value} {est_cost.unit} "
                                "(incomplete: at least one leaf lacks listed price)"
                            )
                        else:
                            basis_desc = (
                                f"step_yield={default_step_yield:g}, stoich_excess={stoich_excess:g}"
                                if (default_step_yield < 1.0 or stoich_excess != 1.0)
                                else "per 1 g target; ignores solvents, yield, labor"
                            )
                            reasons.append(
                                f"estimated leaf-only cost: {est_cost.value} {est_cost.unit} "
                                f"({basis_desc})"
                            )
                        skipped_note = skipped_leaves_from_method(best_route.est_cost_method)
                        if skipped_note:
                            reasons.append(f"leaves excluded from route cost: {skipped_note}")
                    else:
                        status = "unknown"
                        n_steps = None
                        best_price = None
                        est_cost = None
                        cost_complete = False
                        reasons.append(
                            f"not found in stock backends and no complete route to purchasable leaves "
                            f"found within max_steps={max_steps} (unknown does not mean unmakeable)"
                        )
                        if routes:
                            reasons.append(
                                f"found {len(routes)} incomplete candidate disconnection(s) with "
                                "unpurchasable or unverified leaves"
                            )

                used_sources, lics = _collect_provenance_and_licenses(
                    offers, routes, sa_scores, leaf_offer_cache
                )
                entry = _EvaluatedEntry(
                    canonical_smiles=nm.canonical_smiles,
                    inchikey=nm.inchikey,
                    formula=nm.formula,
                    status=status,
                    n_steps=n_steps,
                    best_price=best_price,
                    est_cost=est_cost,
                    est_cost_yield_adjusted=est_cost_yield_adjusted,
                    min_outlay_cost=min_outlay_cost,
                    step_yield_assumed=default_step_yield,
                    cost_complete=cost_complete,
                    lead_time_days=lead_time_days,
                    moq_pack_g=moq_pack_g,
                    precursor_smiles=precursor_smiles,
                    inchikey_match=inchikey_match,
                    matched_offer_inchikeys=matched_offer_iks,
                    stereoisomer_inchikeys_in_stock=stereo_iks,
                    stereoisomer_offers=list(stereo_offers),
                    leaf_offers=rec_leaf_offers,
                    sa_scores=sa_scores,
                    offers=offers,
                    routes=routes,
                    reasons=reasons,
                    evidence=evidence,
                    used_sources=used_sources,
                    licenses=lics,
                )
                memo_by_ik[nm.inchikey] = entry

            mol_obj = Molecule(
                smiles=entry.canonical_smiles,
                inchikey=entry.inchikey,
                name=cand_name,
                formula=entry.formula,
            )
            rec = AvailabilityRecord(
                molecule=mol_obj,
                offers=list(entry.offers),
                routes=list(entry.routes),
                leaf_offers={k: list(v) for k, v in entry.leaf_offers.items()},
                lead_time_days=entry.lead_time_days,
                moq_pack_g=entry.moq_pack_g,
                precursor_smiles=list(entry.precursor_smiles),
                inchikey_match=entry.inchikey_match,
                matched_offer_inchikeys=list(entry.matched_offer_inchikeys),
                stereoisomer_offers=list(entry.stereoisomer_offers),
                sa_scores=entry.sa_scores,
                origin="real",
                provenance=Provenance(
                    source_url=entry.evidence[0] if entry.evidence else None,
                    software=[SoftwareRef(name="buy-or-make", version="0.1.0")],
                    used=list(entry.used_sources),
                ),
                licenses=list(entry.licenses),
            ).sealed()
            records.append(rec)

            results.append(
                RankResult(
                    candidate=raw,
                    name=cand_name,
                    canonical_smiles=entry.canonical_smiles,
                    inchikey=entry.inchikey,
                    status=entry.status,  # type: ignore[arg-type]
                    n_steps=entry.n_steps,
                    best_price=entry.best_price,
                    est_cost=entry.est_cost,
                    est_cost_yield_adjusted=entry.est_cost_yield_adjusted,
                    min_outlay_cost=entry.min_outlay_cost,
                    step_yield_assumed=entry.step_yield_assumed,
                    cost_complete=entry.cost_complete,
                    lead_time_days=entry.lead_time_days,
                    moq_pack_g=entry.moq_pack_g,
                    precursor_smiles=list(entry.precursor_smiles),
                    inchikey_match=entry.inchikey_match,
                    matched_offer_inchikeys=list(entry.matched_offer_inchikeys),
                    stereoisomer_inchikeys_in_stock=list(entry.stereoisomer_inchikeys_in_stock),
                    sa_scores=entry.sa_scores,
                    reasons=list(entry.reasons),
                    evidence=list(entry.evidence),
                )
            )
    finally:
        cache.close()
        if local_stock is not None:
            local_stock.close()
        if pubchem is not None:
            pubchem.close()
        lookup_backend.close()

    results.sort(key=_rank_sort_key)
    return results, records


def cost_aware_rank(
    candidates: list[str],
    max_steps: int = 3,
    backends: tuple[str, ...] = ("local_stock", "pubchem"),
    route_backend: str = "auto",
    offline: bool = False,
    **kwargs,
) -> list[RankResult]:
    """Rank candidate SMILES: buyable (cheaper known price first) < makeable (fewer steps, lower cost) < unknown < invalid."""
    ranked, _ = evaluate_candidates(
        candidates,
        max_steps=max_steps,
        backends=backends,
        route_backend=route_backend,
        offline=offline,
        **kwargs,
    )
    return ranked


def filter_decisions_from_ranked(
    ranked: Iterable[RankResult],
    *,
    budget: Optional[float] = None,
    budget_currency: str = "USD",
    require_complete_cost: bool = False,
    allow_raw_pack_budget: bool = False,
    max_outlay_budget: Optional[float] = None,
    enforce_moq_outlay: bool = False,
    use_yield_adjusted_cost: bool = False,
) -> FilterSummary:
    """Apply unit/currency-aware budget, yield-basis, and discrete MOQ outlay gating to ranked results."""
    kept: list[RankResult] = []
    over_budget_list: list[RankResult] = []
    rejected: list[RankResult] = []
    decisions: list[ProposalDecision] = []

    per_g_unit = f"{budget_currency}/g"
    buyable_allowed_units = {per_g_unit, budget_currency} if allow_raw_pack_budget else {per_g_unit}
    makeable_allowed_units = {budget_currency, per_g_unit}
    eff_outlay_limit = (
        max_outlay_budget
        if max_outlay_budget is not None
        else (budget if enforce_moq_outlay else None)
    )

    for res in ranked:
        reasons = list(res.reasons)
        is_actionable = res.status in {"buyable", "makeable"}
        is_over_budget = False
        price_status: PriceStatusLiteral = "not_applicable"
        keep = is_actionable

        if is_actionable and budget is not None:
            if res.status == "buyable":
                relevant_cost = res.best_price
                allowed_units = buyable_allowed_units
            else:
                relevant_cost = (
                    res.est_cost_yield_adjusted
                    if (use_yield_adjusted_cost and res.est_cost_yield_adjusted is not None)
                    else res.est_cost
                )
                allowed_units = makeable_allowed_units

            if relevant_cost is None or relevant_cost.value is None:
                price_status = "missing"
                reasons.append(
                    f"missing_price: listed price is null; cannot verify <= {budget} {per_g_unit}"
                )
                if require_complete_cost:
                    keep = False
                    reasons.append("rejected: require_complete_cost=True and price is missing")
            elif relevant_cost.unit not in allowed_units:
                price_status = "incompatible_unit"
                keep = False
                reasons.append(
                    f"incompatible_cost_unit: offer unit '{relevant_cost.unit}' is not comparable "
                    f"to '{per_g_unit}' budget"
                )
            elif relevant_cost.value > budget:
                is_over_budget = True
                keep = False
                price_status = (
                    "incomplete"
                    if (res.status == "makeable" and not res.cost_complete)
                    else "verified"
                )
                reasons.append(
                    f"over_budget: cost {relevant_cost.value} {relevant_cost.unit} > budget {budget} {per_g_unit}"
                )
            else:
                if res.status == "makeable" and not res.cost_complete:
                    price_status = "incomplete"
                    reasons.append(
                        f"incomplete_cost: partial route cost {relevant_cost.value} {relevant_cost.unit} "
                        f"<= budget {budget} {per_g_unit}, but unpriced leaves remain"
                    )
                    if require_complete_cost:
                        keep = False
                        reasons.append(
                            "rejected: require_complete_cost=True and route cost is incomplete"
                        )
                elif (
                    res.status == "makeable"
                    and require_complete_cost
                    and res.step_yield_assumed >= 1.0
                    and not use_yield_adjusted_cost
                ):
                    price_status = "leaf_lower_bound"
                    reasons.append(
                        "leaf_lower_bound: route cost assumes 100% step yield and fractional pack mass"
                    )
                else:
                    price_status = "verified"

        if keep and is_actionable and eff_outlay_limit is not None:
            outlay = res.min_outlay_cost
            if outlay is None or outlay.value is None:
                if require_complete_cost:
                    keep = False
                    price_status = "missing"
                    reasons.append(
                        f"missing_moq_outlay: discrete pack cash outlay unknown; "
                        f"cannot verify <= {eff_outlay_limit} {budget_currency}"
                    )
            elif outlay.unit != budget_currency:
                keep = False
                price_status = "incompatible_unit"
                reasons.append(
                    f"incompatible_outlay_unit: '{outlay.unit}' != '{budget_currency}'"
                )
            elif outlay.value > eff_outlay_limit:
                keep = False
                is_over_budget = True
                reasons.append(
                    f"over_budget_moq_outlay: discrete pack cash outlay {outlay.value} {outlay.unit} "
                    f"> outlay budget {eff_outlay_limit} {budget_currency}"
                )

        if keep:
            kept.append(res)
        elif is_over_budget:
            over_budget_list.append(res)
            rejected.append(res)
        else:
            rejected.append(res)

        decisions.append(
            ProposalDecision(
                candidate=res.candidate,
                keep=keep,
                over_budget=is_over_budget,
                price_status=price_status,
                status=res.status,
                reasons=reasons,
                result=res,
            )
        )

    return FilterSummary(
        budget=budget,
        budget_currency=budget_currency,
        require_complete_cost=require_complete_cost,
        max_outlay_budget=eff_outlay_limit,
        use_yield_adjusted_cost=use_yield_adjusted_cost,
        kept=kept,
        over_budget=over_budget_list,
        rejected=rejected,
        decisions=decisions,
    )


def filter_proposals(
    proposals: list[str],
    max_steps: int = 3,
    budget: Optional[float] = None,
    *,
    budget_currency: str = "USD",
    require_complete_cost: bool = False,
    allow_raw_pack_budget: bool = False,
    max_outlay_budget: Optional[float] = None,
    enforce_moq_outlay: bool = False,
    use_yield_adjusted_cost: bool = False,
    **kwargs,
) -> FilterSummary:
    """Filter co-scientist proposals: keep buyable/makeable within budget, enforcing unit & completeness policy."""
    kwargs.setdefault("currency", budget_currency)
    ranked = cost_aware_rank(proposals, max_steps=max_steps, **kwargs)
    return filter_decisions_from_ranked(
        ranked,
        budget=budget,
        budget_currency=budget_currency,
        require_complete_cost=require_complete_cost,
        allow_raw_pack_budget=allow_raw_pack_budget,
        max_outlay_budget=max_outlay_budget,
        enforce_moq_outlay=enforce_moq_outlay,
        use_yield_adjusted_cost=use_yield_adjusted_cost,
    )
