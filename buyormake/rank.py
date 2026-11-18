"""Cost-aware ranking and proposal filtering for candidate molecules."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

import httpx

from .cache import DEFAULT_CACHE_PATH, ResponseCache
from .canonical import normalize_smiles
from .cost import cheapest_offer_for_target
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
    ProposalDecision,
    Provenance,
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


def _rank_sort_key(res: RankResult) -> tuple:
    if res.status == "buyable":
        has_price = res.best_price is not None and res.best_price.value is not None
        price_val = float(res.best_price.value) if has_price else 0.0
        return (0, 0 if has_price else 1, price_val)
    if res.status == "makeable":
        has_cost = res.est_cost is not None and res.est_cost.value is not None
        if has_cost and res.cost_complete:
            cost_tier = 0
        elif has_cost:
            cost_tier = 1
        else:
            cost_tier = 2
        cost_val = float(res.est_cost.value) if has_cost else 0.0
        return (1, res.n_steps if res.n_steps is not None else 999, cost_tier, cost_val)
    if res.status == "unknown":
        return (2, 0, 0, 0.0)
    return (3, 0, 0, 0.0)


def _unique_urls(offers: Iterable[Offer]) -> list[str]:
    seen: set[str] = set()
    urls: list[str] = []
    for off in offers:
        if off.source_url and off.source_url not in seen:
            seen.add(off.source_url)
            urls.append(off.source_url)
    return urls


def evaluate_candidates(
    candidates: list[str],
    *,
    max_steps: int = 3,
    backends: tuple[str, ...] = ("local_stock", "pubchem"),
    route_backend: str = "auto",
    offline: bool = False,
    cache_path: Path | str = DEFAULT_CACHE_PATH,
    stock_paths: Optional[list[Path | str]] = None,
    reaction_db_path: Path | str = DEFAULT_REACTION_DB,
    scscore_weights: Path | str = DEFAULT_SCSCORE_WEIGHTS,
    aizynth_config: Optional[Path | str] = None,
    include_scores: bool = True,
    http_client: Optional[httpx.Client] = None,
) -> tuple[list[RankResult], list[AvailabilityRecord]]:
    """Evaluate candidates and return (sorted RankResults, AvailabilityRecords in input order)."""
    cache = ResponseCache(cache_path)
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

    try:
        for raw in candidates:
            nm = normalize_smiles(raw)
            if nm is None:
                results.append(
                    RankResult(
                        candidate=raw,
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

            offers = lookup_offers(nm.inchikey, nm.canonical_smiles)
            sa_scores = (
                compute_scores(nm.canonical_smiles, scscore_weights=scscore_weights)
                if include_scores
                else SAScores()
            )
            routes: list[Route] = []
            reasons: list[str] = []
            evidence = _unique_urls(offers)

            if offers:
                status = "buyable"
                best_price = cheapest_offer_for_target(offers)
                est_cost = None
                cost_complete = best_price is not None
                n_steps = 0
                vendors_preview = ", ".join(dict.fromkeys(o.vendor for o in offers[:3]))
                suffix = f" (+{len(offers) - 3} more offers)" if len(offers) > 3 else ""
                reasons.append(f"listed by {len(offers)} catalog offer(s): {vendors_preview}{suffix}")
                if best_price is None:
                    reasons.append("listed price unknown (source does not publish open catalog prices)")
                else:
                    reasons.append(f"cheapest listed price: {best_price.raw_unit_str or best_price.value}")
            else:
                routes = find_routes(
                    nm.canonical_smiles,
                    nm.inchikey,
                    max_steps=max_steps,
                    route_backend=route_backend,
                    offer_lookup=lookup_offers,
                    lookup_backend=lookup_backend,
                    template_backend=template_backend,
                    aizynth_backend=aizynth_backend,
                    offer_cache=offer_cache,
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

            used_sources = sorted({o.source_id for o in offers})
            lics = sorted({SOURCE_BY_ID[s].license for s in used_sources if s in SOURCE_BY_ID})
            if sa_scores.sascore is not None:
                lics.append("BSD-3-Clause")
            if sa_scores.scscore is not None:
                lics.append("MIT")
            lics = sorted(set(lics))

            mol_obj = Molecule(
                smiles=nm.canonical_smiles,
                inchikey=nm.inchikey,
                formula=nm.formula,
            )
            rec = AvailabilityRecord(
                molecule=mol_obj,
                offers=offers,
                routes=routes,
                sa_scores=sa_scores,
                origin="real",
                provenance=Provenance(
                    source_url=evidence[0] if evidence else None,
                    software=[SoftwareRef(name="buy-or-make", version="0.1.0")],
                    used=used_sources,
                ),
                licenses=lics,
            ).sealed()
            records.append(rec)

            results.append(
                RankResult(
                    candidate=raw,
                    canonical_smiles=nm.canonical_smiles,
                    inchikey=nm.inchikey,
                    status=status,
                    n_steps=n_steps,
                    best_price=best_price,
                    est_cost=est_cost,
                    cost_complete=cost_complete,
                    sa_scores=sa_scores,
                    reasons=reasons,
                    evidence=evidence,
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


def filter_proposals(
    proposals: list[str],
    max_steps: int = 3,
    budget: Optional[float] = None,
    **kwargs,
) -> FilterSummary:
    """Filter co-scientist proposals: keep buyable/makeable, flag over-budget, return reasons."""
    ranked = cost_aware_rank(proposals, max_steps=max_steps, **kwargs)
    kept: list[RankResult] = []
    over_budget_list: list[RankResult] = []
    rejected: list[RankResult] = []
    decisions: list[ProposalDecision] = []

    for res in ranked:
        reasons = list(res.reasons)
        is_actionable = res.status in {"buyable", "makeable"}
        is_over_budget = False

        if is_actionable and budget is not None:
            relevant_cost = res.best_price if res.status == "buyable" else res.est_cost
            if relevant_cost is not None and relevant_cost.value is not None:
                if relevant_cost.value > budget:
                    is_over_budget = True
                    reasons.append(
                        f"over_budget: cost {relevant_cost.value} {relevant_cost.unit} > budget {budget}"
                    )

        keep = is_actionable and not is_over_budget
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
                status=res.status,
                reasons=reasons,
                result=res,
            )
        )

    return FilterSummary(
        budget=budget,
        kept=kept,
        over_budget=over_budget_list,
        rejected=rejected,
        decisions=decisions,
    )
