"""Cost-aware ranking and proposal filtering for candidate molecules."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import httpx

from .cache import DEFAULT_CACHE_PATH, DEFAULT_ONLINE_TTL_SECONDS, ResponseCache
from .canonical import normalize_smiles
from .cost import cheapest_offer_for_target, skipped_leaves_from_method
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
}


def _rank_sort_key(res: RankResult) -> tuple[int, int, int, str, float]:
    if res.status == "buyable":
        if res.best_price is not None and res.best_price.value is not None:
            price_val = float(res.best_price.value)
            unit = res.best_price.unit
            if unit == "USD/g":
                return (0, 0, 0, unit, price_val)
            # Group other currencies (EUR/g, JPY/g) and raw per-pack prices (USD) separately by unit
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
            mapped_src = _STEP_SOURCE_TO_SOURCE_ID.get(st.source)
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
    cost_complete: bool
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
    cache = ResponseCache(cache_path, max_age_seconds=effective_ttl, refresh=refresh_cache)
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

    offer_cache: dict[str, list[Offer]] = {}

    def lookup_offers(inchikey: str, canonical_smiles: str) -> list[Offer]:
        if inchikey in offer_cache:
            return offer_cache[inchikey]
        found: list[Offer] = []
        if local_stock is not None:
            found.extend(local_stock.lookup(inchikey))
        if pubchem is not None:
            found.extend(pubchem.lookup(inchikey, canonical_smiles))
        offer_cache[inchikey] = found
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
                        cost_complete=False,
                        sa_scores=None,
                        reasons=["invalid SMILES; RDKit could not parse structure"],
                        evidence=[],
                    )
                )
                continue

            entry = memo_by_ik.get(nm.inchikey)
            if entry is None:
                offers = lookup_offers(nm.inchikey, nm.canonical_smiles)
                desalt_reason: Optional[str] = None
                parent_nm = None
                if not offers and desalt_fallback and (
                    "." in nm.canonical_smiles
                    or "+" in nm.canonical_smiles
                    or "-" in nm.canonical_smiles
                ):
                    cand_parent = normalize_smiles(raw, desalt=True)
                    if cand_parent is not None and cand_parent.inchikey != nm.inchikey:
                        parent_nm = cand_parent
                        parent_offers = lookup_offers(
                            parent_nm.inchikey, parent_nm.canonical_smiles
                        )
                        if parent_offers:
                            offers = parent_offers
                            desalt_reason = (
                                f"matched desalted parent structure {parent_nm.canonical_smiles} "
                                f"({parent_nm.inchikey}) after exact salt/mixture form had 0 offers"
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
                reasons: list[str] = []
                if desalt_reason is not None:
                    reasons.append(desalt_reason)
                evidence = _unique_urls(offers)

                if offers:
                    status = "buyable"
                    best_price = cheapest_offer_for_target(offers, currency=currency)
                    est_cost = None
                    # cost_complete describes a route cost estimate; buyables have none (see best_price)
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
                        offer_lookup=lookup_offers,
                        lookup_backend=lookup_backend,
                        template_backend=template_backend,
                        aizynth_backend=aizynth_backend,
                        offer_cache=offer_cache,
                        target_mol=route_target.mol,
                        currency=currency,
                    )
                    solved = [r for r in routes if r.leaves and all(r.leaves_purchasable)]
                    if solved:
                        best_route = solved[0]
                        status = "makeable"
                        n_steps = best_route.n_steps
                        best_price = None
                        est_cost = best_route.est_cost
                        cost_complete = best_route.est_cost_complete
                        leaf_urls: list[str] = []
                        for leaf_smi in best_route.leaves:
                            lnm = normalize_smiles(leaf_smi)
                            if lnm is not None:
                                leaf_urls.extend(_unique_urls(offer_cache.get(lnm.inchikey, [])))
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
                            reasons.append(
                                f"estimated leaf-only cost: {est_cost.value} {est_cost.unit} "
                                "(per 1 g target; ignores solvents, yield, labor)"
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
                    offers, routes, sa_scores, offer_cache
                )
                entry = _EvaluatedEntry(
                    canonical_smiles=nm.canonical_smiles,
                    inchikey=nm.inchikey,
                    formula=nm.formula,
                    status=status,
                    n_steps=n_steps,
                    best_price=best_price,
                    est_cost=est_cost,
                    cost_complete=cost_complete,
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
                    cost_complete=entry.cost_complete,
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
) -> FilterSummary:
    """Apply unit/currency-aware budget and price-completeness gating to ranked results."""
    kept: list[RankResult] = []
    over_budget_list: list[RankResult] = []
    rejected: list[RankResult] = []
    decisions: list[ProposalDecision] = []

    per_g_unit = f"{budget_currency}/g"
    buyable_allowed_units = {per_g_unit, budget_currency} if allow_raw_pack_budget else {per_g_unit}
    makeable_allowed_units = {budget_currency, per_g_unit}

    for res in ranked:
        reasons = list(res.reasons)
        is_actionable = res.status in {"buyable", "makeable"}
        is_over_budget = False
        price_status: PriceStatusLiteral = "not_applicable"
        keep = is_actionable

        if is_actionable and budget is not None:
            relevant_cost = res.best_price if res.status == "buyable" else res.est_cost
            allowed_units = (
                buyable_allowed_units if res.status == "buyable" else makeable_allowed_units
            )

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
                else:
                    price_status = "verified"

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
    )
