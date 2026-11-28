"""Honest route cost estimator from real listed catalog prices only, with step-yield and MOQ bounds."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Optional

from .canonical import normalize_smiles
from .schema import Offer, Quantity, RouteStep

MASS_TO_GRAMS = {
    "ug": 1e-6,
    "\u00b5g": 1e-6,
    "mg": 1e-3,
    "g": 1.0,
    "kg": 1e3,
}

DEFAULT_REALISTIC_STEP_YIELD = 0.75
DEFAULT_REALISTIC_STOICH_EXCESS = 1.05

EST_COST_METHOD = (
    "sum over leaves of (cheapest listed USD/g * 1 equiv mass per 1 g target, "
    "i.e., count * MW_leaf / MW_target); ignores solvents, reagents, catalysts, "
    "yield, labor, and purification"
)


def purchasable_offers(offers: Iterable[Offer]) -> list[Offer]:
    """Return offers not explicitly flagged out of stock (`in_stock is not False`).

    Offers with `in_stock=True` or `in_stock=None` (stock flag unpublished by free
    catalogs) remain purchasable; only explicit `in_stock=False` is excluded.
    """
    return [off for off in offers if off.in_stock is not False]


def offer_pack_grams(offer: Offer) -> Optional[float]:
    """Return pack mass in grams if `offer.pack_size` carries a recognized mass unit."""
    if offer.pack_size is None or offer.pack_size.value is None or offer.pack_size.value <= 0:
        return None
    factor = MASS_TO_GRAMS.get(offer.pack_size.unit)
    if factor is None:
        return None
    grams = offer.pack_size.value * factor
    return grams if grams > 0 else None


def offer_price_per_gram(
    offer: Offer,
    currency: str = "USD",
    *,
    max_pack_grams: Optional[float] = None,
) -> Optional[float]:
    """Return USD/g (or matching currency/g) if offer has both listed price and mass pack_size."""
    if offer.price is None or offer.currency != currency:
        return None
    grams = offer_pack_grams(offer)
    if grams is None:
        return None
    if max_pack_grams is not None and grams > max_pack_grams + 1e-12:
        return None
    return offer.price / grams


def _molar_mass_ratio(
    offer_ik: str,
    target_ik: Optional[str],
    target_mw: Optional[float],
    offer_mol_weights: Optional[Mapping[str, float]],
) -> tuple[float, Optional[float]]:
    """Return `(mw_offer / mw_target, mw_offer)` when `offer_ik` is a salt/parent match to `target_ik`."""
    if (
        target_ik is None
        or offer_ik == target_ik
        or target_mw is None
        or target_mw <= 0
        or not offer_mol_weights
    ):
        return 1.0, None
    mw_off = offer_mol_weights.get(offer_ik)
    if mw_off is None or mw_off <= 0 or abs(mw_off - target_mw) <= 1e-3:
        return 1.0, None
    return float(mw_off) / float(target_mw), float(mw_off)


def _pack_key(off: Offer) -> Optional[tuple]:
    if off.pack_size is None or off.pack_size.value is None:
        return None
    return (off.pack_size.value, off.pack_size.unit)


def _cheapest_in_currency(
    offers: list[Offer],
    currency: str,
    *,
    target_mol_weight: Optional[float] = None,
    target_inchikey: Optional[str] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
    max_pack_grams: Optional[float] = None,
) -> Optional[Quantity]:
    priced_per_g: list[tuple[float, Offer, float, Optional[float]]] = []
    priced_raw: list[tuple[float, Offer]] = []
    for off in offers:
        if off.price is None or off.currency != currency:
            continue
        ppg = offer_price_per_gram(off, currency=currency, max_pack_grams=max_pack_grams)
        if ppg is not None:
            mw_ratio, mw_off = _molar_mass_ratio(
                off.inchikey, target_inchikey, target_mol_weight, offer_mol_weights
            )
            priced_per_g.append((ppg * mw_ratio, off, mw_ratio, mw_off))
        elif max_pack_grams is None:
            priced_raw.append((off.price, off))

    if priced_per_g:
        adj_ppg, off, mw_ratio, mw_off = min(priced_per_g, key=lambda t: t[0])
        pack_str = off.pack_size.raw_unit_str if off.pack_size else ""
        base_str = f"{off.price} {off.currency} / {pack_str}".strip()
        if mw_off is not None and target_mol_weight is not None:
            base_str = (
                f"{base_str} (molar-adjusted x{mw_ratio:.4f} for "
                f"MW {mw_off:.2f}/{target_mol_weight:.2f})"
            )
        return Quantity(
            value=round(adj_ppg, 4),
            unit=f"{currency}/g",
            raw_unit_str=base_str,
        )
    if priced_raw:
        pack_keys = {_pack_key(off) for _, off in priced_raw}
        comparable = len(priced_raw) == 1 or (len(pack_keys) == 1 and None not in pack_keys)
        if comparable:
            val, off = min(priced_raw, key=lambda t: t[0])
            pack_str = (
                off.pack_size.raw_unit_str
                if (off.pack_size and off.pack_size.raw_unit_str)
                else "pack"
            )
            return Quantity(
                value=round(val, 4),
                unit=currency,
                raw_unit_str=f"{val} {currency} / {pack_str}",
            )
    return None


def cheapest_offer_for_target(
    offers: list[Offer],
    currency: str = "USD",
    *,
    target_mol_weight: Optional[float] = None,
    target_inchikey: Optional[str] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
    max_pack_grams: Optional[float] = None,
) -> Optional[Quantity]:
    """Return cheapest listed price among offers, normalized to currency/g when mass pack_size is known.

    When an offer has a salt/parent InChIKey discrepancy relative to `target_inchikey` and molecular
    weights are supplied via `target_mol_weight` + `offer_mol_weights`, scales `currency/g` by
    `MW_offer / MW_target` to preserve equimolar active-species stoichiometry.
    """
    found = _cheapest_in_currency(
        offers,
        currency,
        target_mol_weight=target_mol_weight,
        target_inchikey=target_inchikey,
        offer_mol_weights=offer_mol_weights,
        max_pack_grams=max_pack_grams,
    )
    if found is not None:
        return found
    other_currencies = {
        off.currency for off in offers if off.price is not None and off.currency is not None
    }
    if len(other_currencies) == 1:
        (only,) = other_currencies
        if only != currency:
            return _cheapest_in_currency(
                offers,
                only,
                target_mol_weight=target_mol_weight,
                target_inchikey=target_inchikey,
                offer_mol_weights=offer_mol_weights,
                max_pack_grams=max_pack_grams,
            )
    return None


def min_outlay_for_offers(
    offers: list[Offer],
    currency: str = "USD",
    *,
    required_grams: float = 1.0,
    target_mol_weight: Optional[float] = None,
    target_inchikey: Optional[str] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
) -> tuple[Optional[Quantity], Optional[Offer], Optional[float]]:
    """Return `(min_outlay_quantity, selected_offer, moq_pack_g)` for procuring `required_grams`.

    Minimizes discrete pack cash outlay `ceil(required_offer_grams / pack_grams) * offer.price`.
    """
    if required_grams <= 0 or not math.isfinite(required_grams):
        raise ValueError("required_grams must be a finite number > 0")

    best_priced: Optional[tuple[float, float, Offer, float]] = None
    best_unpriced_pack: Optional[tuple[float, Offer]] = None

    for off in offers:
        pack_g = offer_pack_grams(off)
        if pack_g is None:
            continue
        mw_ratio, _ = _molar_mass_ratio(
            off.inchikey, target_inchikey, target_mol_weight, offer_mol_weights
        )
        eff_req_g = required_grams * mw_ratio
        if off.price is not None and off.currency == currency:
            n_packs = max(1, int(math.ceil(eff_req_g / pack_g - 1e-12)))
            outlay = n_packs * off.price
            ppg = off.price / pack_g
            cand = (outlay, ppg, off, pack_g)
            if best_priced is None or (cand[0], cand[1]) < (best_priced[0], best_priced[1]):
                best_priced = cand
        else:
            if best_unpriced_pack is None or pack_g < best_unpriced_pack[0]:
                best_unpriced_pack = (pack_g, off)

    if best_priced is not None:
        outlay, _, sel_off, pack_g = best_priced
        pack_str = sel_off.pack_size.raw_unit_str if sel_off.pack_size else f"{pack_g:g} g"
        return (
            Quantity(
                value=round(outlay, 4),
                unit=currency,
                raw_unit_str=f"{round(outlay, 4)} {currency} discrete pack MOQ outlay ({pack_str})",
            ),
            sel_off,
            round(pack_g, 6),
        )
    if best_unpriced_pack is not None:
        pack_g, sel_off = best_unpriced_pack
        return None, sel_off, round(pack_g, 6)
    return None, (offers[0] if offers else None), None


def compute_leaf_depths(
    target_smiles: str,
    leaves: list[str],
    steps: Optional[list[RouteStep]] = None,
) -> list[tuple[str, int]]:
    """Compute exact reaction-tree step depth `d_i >= 1` for each leaf occurrence in `leaves`.

    Distinguishes deep linear sequences (early leaves at `d = n_steps`) from convergent
    sequences of equal step count (parallel branches at shallower `d`).
    """
    if not leaves:
        return []
    if not steps:
        return [(leaf, 1) for leaf in leaves]

    def _can(smi: str) -> str:
        nm = normalize_smiles(smi)
        return nm.canonical_smiles if nm is not None else smi.strip()

    root_can = _can(target_smiles)
    steps_by_prod: dict[str, list[RouteStep]] = {}
    first_step_by_prod: dict[str, RouteStep] = {}
    for st in steps:
        p_can = _can(st.product)
        steps_by_prod.setdefault(p_can, []).append(st)
        first_step_by_prod.setdefault(p_can, st)

    walked_by_leaf: dict[str, list[int]] = {}

    def _walk(mol_can: str, depth: int, ancestors: frozenset[str]) -> None:
        if mol_can not in ancestors:
            queue = steps_by_prod.get(mol_can)
            st: Optional[RouteStep] = None
            if queue:
                st = queue.pop(0)
            elif mol_can in first_step_by_prod:
                st = first_step_by_prod[mol_can]
            if st is not None and st.reactants:
                next_anc = ancestors | {mol_can}
                for r_smi in st.reactants:
                    _walk(_can(r_smi), depth + 1, next_anc)
                return
        if depth > 0:
            walked_by_leaf.setdefault(mol_can, []).append(depth)

    _walk(root_can, 0, frozenset())

    fallback_depth = max(1, len(steps))
    result: list[tuple[str, int]] = []
    for leaf_smi in leaves:
        l_can = _can(leaf_smi)
        depths = walked_by_leaf.get(l_can)
        if depths:
            d = depths.pop(0)
        else:
            d = fallback_depth
        result.append((leaf_smi, max(1, int(d))))
    return result


@dataclass(frozen=True)
class RouteCostBreakdown:
    """Detailed stoichiometric, yield-compounded, and discrete pack MOQ route cost bounds."""

    est_cost: Optional[Quantity]
    est_cost_complete: bool
    lower_bound_cost: Optional[Quantity]
    yield_adjusted_cost: Optional[Quantity]
    min_outlay_cost: Optional[Quantity]
    min_outlay_complete: bool
    moq_pack_g: Optional[float]
    leaf_depths: dict[str, list[int]] = field(default_factory=dict)
    skipped: list[str] = field(default_factory=list)


def _validate_cost_params(
    default_step_yield: float,
    stoich_excess: float,
    target_scale_g: float,
    max_pack_multiple: Optional[float],
) -> None:
    if (
        isinstance(default_step_yield, bool)
        or not math.isfinite(default_step_yield)
        or not (0.0 < default_step_yield <= 1.0)
    ):
        raise ValueError("default_step_yield must be in (0, 1]")
    if isinstance(stoich_excess, bool) or not math.isfinite(stoich_excess) or stoich_excess <= 0.0:
        raise ValueError("stoich_excess must be > 0")
    if isinstance(target_scale_g, bool) or not math.isfinite(target_scale_g) or target_scale_g <= 0.0:
        raise ValueError("target_scale_g must be > 0")
    if max_pack_multiple is not None and (
        isinstance(max_pack_multiple, bool)
        or not math.isfinite(max_pack_multiple)
        or max_pack_multiple <= 0.0
    ):
        raise ValueError("max_pack_multiple must be > 0 when set")


def estimate_route_cost_bounds(
    target_smiles: str,
    leaves: list[str],
    leaf_offers: dict[str, list[Offer]],
    currency: str = "USD",
    *,
    steps: Optional[list[RouteStep]] = None,
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
    realistic_step_yield: float = DEFAULT_REALISTIC_STEP_YIELD,
    realistic_stoich_excess: float = DEFAULT_REALISTIC_STOICH_EXCESS,
    target_scale_g: float = 1.0,
    max_pack_multiple: Optional[float] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
) -> RouteCostBreakdown:
    """Compute primary, 100%-yield lower-bound, yield-compounded, and discrete pack MOQ costs."""
    _validate_cost_params(default_step_yield, stoich_excess, target_scale_g, max_pack_multiple)
    _validate_cost_params(realistic_step_yield, realistic_stoich_excess, target_scale_g, None)

    if not leaves:
        return RouteCostBreakdown(
            est_cost=None,
            est_cost_complete=False,
            lower_bound_cost=None,
            yield_adjusted_cost=None,
            min_outlay_cost=None,
            min_outlay_complete=False,
            moq_pack_g=None,
        )

    target_mol = normalize_smiles(target_smiles)
    if target_mol is None or target_mol.mol_weight <= 0:
        return RouteCostBreakdown(
            est_cost=None,
            est_cost_complete=False,
            lower_bound_cost=None,
            yield_adjusted_cost=None,
            min_outlay_cost=None,
            min_outlay_complete=False,
            moq_pack_g=None,
            skipped=[f"{target_smiles} (target unparseable; cost not estimated)"],
        )

    leaf_depth_pairs = compute_leaf_depths(target_smiles, leaves, steps)
    depths_by_leaf: dict[str, list[int]] = {}
    for leaf_smi, depth in leaf_depth_pairs:
        depths_by_leaf.setdefault(leaf_smi, []).append(depth)

    eff_real_yield = default_step_yield if default_step_yield < 1.0 else realistic_step_yield
    eff_real_excess = stoich_excess if stoich_excess > 1.0 else realistic_stoich_excess

    total_primary = 0.0
    total_lower = 0.0
    total_yield_adj = 0.0
    total_outlay = 0.0

    priced_leaves = 0
    outlay_leaves = 0
    total_unique_leaves = len(depths_by_leaf)
    skipped: list[str] = []
    leaf_moq_packs: list[float] = []

    for leaf_smi, occ_depths in depths_by_leaf.items():
        leaf_mol = normalize_smiles(leaf_smi)
        if leaf_mol is None:
            skipped.append(f"{leaf_smi} (unparseable SMILES)")
            continue
        raw_offers = leaf_offers.get(leaf_mol.inchikey, [])
        offers = purchasable_offers(raw_offers)

        # Molar equivalents per mole of target under (1) primary params, (2) ideal 100% lower bound, (3) realistic yield
        primary_equiv = sum(stoich_excess / (default_step_yield**d) for d in occ_depths)
        lower_equiv = float(len(occ_depths))
        real_equiv = sum(eff_real_excess / (eff_real_yield**d) for d in occ_depths)

        mw_frac = leaf_mol.mol_weight / target_mol.mol_weight
        primary_mass_per_g = primary_equiv * mw_frac
        lower_mass_per_g = lower_equiv * mw_frac
        real_mass_per_g = real_equiv * mw_frac
        req_grams_at_scale = target_scale_g * primary_mass_per_g

        max_pack_g = (
            max(target_scale_g, req_grams_at_scale) * max_pack_multiple
            if max_pack_multiple is not None
            else None
        )

        ppg_candidates: list[float] = []
        for off in offers:
            ppg = offer_price_per_gram(off, currency=currency, max_pack_grams=max_pack_g)
            if ppg is not None:
                mw_ratio, _ = _molar_mass_ratio(
                    off.inchikey, leaf_mol.inchikey, leaf_mol.mol_weight, offer_mol_weights
                )
                ppg_candidates.append(ppg * mw_ratio)

        outlay_q, _, moq_g = min_outlay_for_offers(
            offers,
            currency=currency,
            required_grams=max(req_grams_at_scale, 1e-9),
            target_mol_weight=leaf_mol.mol_weight,
            target_inchikey=leaf_mol.inchikey,
            offer_mol_weights=offer_mol_weights,
        )
        if moq_g is not None:
            leaf_moq_packs.append(moq_g)
        if outlay_q is not None and outlay_q.value is not None:
            total_outlay += outlay_q.value
            outlay_leaves += 1

        if not ppg_candidates:
            if raw_offers and not offers:
                skipped.append(f"{leaf_smi} (all listed offers flagged out of stock)")
            elif max_pack_g is not None and any(
                offer_price_per_gram(off, currency=currency) is not None for off in offers
            ):
                skipped.append(
                    f"{leaf_smi} (no {currency}/g pack <= {round(max_pack_g, 4)} g MOQ cap)"
                )
            else:
                skipped.append(f"{leaf_smi} (no listed {currency}/g price)")
            continue

        best_ppg = min(ppg_candidates)
        total_primary += best_ppg * primary_mass_per_g
        total_lower += best_ppg * lower_mass_per_g
        total_yield_adj += best_ppg * real_mass_per_g
        priced_leaves += 1

    moq_bottleneck_g = round(max(leaf_moq_packs), 6) if (len(leaf_moq_packs) == total_unique_leaves and leaf_moq_packs) else None

    if priced_leaves == 0:
        return RouteCostBreakdown(
            est_cost=None,
            est_cost_complete=False,
            lower_bound_cost=None,
            yield_adjusted_cost=None,
            min_outlay_cost=None,
            min_outlay_complete=False,
            moq_pack_g=moq_bottleneck_g,
            leaf_depths=depths_by_leaf,
            skipped=skipped,
        )

    complete = priced_leaves == total_unique_leaves
    outlay_complete = outlay_leaves == total_unique_leaves
    per_g_unit = f"{currency}/g"
    scale_note = (
        f"per 1 g target; step_yield={default_step_yield:g}, stoich_excess={stoich_excess:g}"
        if (default_step_yield < 1.0 or stoich_excess != 1.0)
        else "per 1 g target"
    )
    suffix = scale_note if complete else f"partial lower bound; {scale_note}"

    est_q = Quantity(
        value=round(total_primary, 4),
        unit=per_g_unit,
        raw_unit_str=f"{round(total_primary, 4)} {per_g_unit} ({suffix})",
    )
    lower_q = Quantity(
        value=round(total_lower, 4),
        unit=per_g_unit,
        raw_unit_str=f"{round(total_lower, 4)} {per_g_unit} (100%-yield lower bound; per 1 g target)",
    )
    yield_q = Quantity(
        value=round(total_yield_adj, 4),
        unit=per_g_unit,
        raw_unit_str=(
            f"{round(total_yield_adj, 4)} {per_g_unit} "
            f"(depth-compounded yield={eff_real_yield:g}, excess={eff_real_excess:g}; per 1 g target)"
        ),
    )
    outlay_q_total = (
        Quantity(
            value=round(total_outlay, 4),
            unit=currency,
            raw_unit_str=(
                f"{round(total_outlay, 4)} {currency} "
                f"(discrete pack MOQ outlay for {target_scale_g:g} g target)"
            ),
        )
        if outlay_leaves > 0
        else None
    )
    return RouteCostBreakdown(
        est_cost=est_q,
        est_cost_complete=complete,
        lower_bound_cost=lower_q,
        yield_adjusted_cost=yield_q,
        min_outlay_cost=outlay_q_total,
        min_outlay_complete=outlay_complete,
        moq_pack_g=moq_bottleneck_g,
        leaf_depths=depths_by_leaf,
        skipped=skipped,
    )


def estimate_route_cost(
    target_smiles: str,
    leaves: list[str],
    leaf_offers: dict[str, list[Offer]],
    currency: str = "USD",
    *,
    steps: Optional[list[RouteStep]] = None,
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
    target_scale_g: float = 1.0,
    max_pack_multiple: Optional[float] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
) -> tuple[Optional[Quantity], bool]:
    """Estimate route cost per 1 g of target from real listed leaf prices only."""
    est, complete, _ = estimate_route_cost_detailed(
        target_smiles,
        leaves,
        leaf_offers,
        currency,
        steps=steps,
        default_step_yield=default_step_yield,
        stoich_excess=stoich_excess,
        target_scale_g=target_scale_g,
        max_pack_multiple=max_pack_multiple,
        offer_mol_weights=offer_mol_weights,
    )
    return est, complete


def estimate_route_cost_detailed(
    target_smiles: str,
    leaves: list[str],
    leaf_offers: dict[str, list[Offer]],
    currency: str = "USD",
    *,
    steps: Optional[list[RouteStep]] = None,
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
    target_scale_g: float = 1.0,
    max_pack_multiple: Optional[float] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
) -> tuple[Optional[Quantity], bool, list[str]]:
    """Like `estimate_route_cost`, also returning one note per leaf excluded from the sum."""
    breakdown = estimate_route_cost_bounds(
        target_smiles,
        leaves,
        leaf_offers,
        currency,
        steps=steps,
        default_step_yield=default_step_yield,
        stoich_excess=stoich_excess,
        target_scale_g=target_scale_g,
        max_pack_multiple=max_pack_multiple,
        offer_mol_weights=offer_mol_weights,
    )
    return breakdown.est_cost, breakdown.est_cost_complete, breakdown.skipped


SKIPPED_LEAVES_PREFIX = "; skipped leaves: "


def route_cost_method(
    skipped: list[str],
    *,
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
) -> str:
    """Return EST_COST_METHOD, annotated with yield/excess basis and any excluded leaves."""
    base = EST_COST_METHOD
    if default_step_yield < 1.0 or stoich_excess != 1.0:
        base = (
            f"sum over leaves of (cheapest listed USD/g * depth-compounded equiv mass per 1 g target, "
            f"step_yield={default_step_yield:g}, stoich_excess={stoich_excess:g})"
        )
    if not skipped:
        return base
    return base + SKIPPED_LEAVES_PREFIX + "; ".join(skipped)


def skipped_leaves_from_method(method: str) -> Optional[str]:
    """Inverse of route_cost_method: return the skipped-leaves note, or None if nothing was skipped."""
    if SKIPPED_LEAVES_PREFIX not in method:
        return None
    return method.split(SKIPPED_LEAVES_PREFIX, 1)[1]
