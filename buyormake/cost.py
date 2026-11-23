"""Honest route cost estimator from real listed catalog prices only."""

from __future__ import annotations

from collections import Counter
from typing import Optional

from .canonical import normalize_smiles
from .schema import Offer, Quantity

MASS_TO_GRAMS = {
    "ug": 1e-6,
    "µg": 1e-6,
    "mg": 1e-3,
    "g": 1.0,
    "kg": 1e3,
}

EST_COST_METHOD = (
    "sum over leaves of (cheapest listed USD/g * 1 equiv mass per 1 g target, "
    "i.e., count * MW_leaf / MW_target); ignores solvents, reagents, catalysts, "
    "yield, labor, and purification"
)


def offer_price_per_gram(offer: Offer, currency: str = "USD") -> Optional[float]:
    """Return USD/g (or matching currency/g) if offer has both listed price and mass pack_size."""
    if offer.price is None or offer.currency != currency:
        return None
    if offer.pack_size is None or offer.pack_size.value is None or offer.pack_size.value <= 0:
        return None
    factor = MASS_TO_GRAMS.get(offer.pack_size.unit)
    if factor is None:
        return None
    grams = offer.pack_size.value * factor
    if grams <= 0:
        return None
    return offer.price / grams


def _pack_key(off: Offer) -> Optional[tuple]:
    if off.pack_size is None or off.pack_size.value is None:
        return None
    return (off.pack_size.value, off.pack_size.unit)


def _cheapest_in_currency(offers: list[Offer], currency: str) -> Optional[Quantity]:
    priced_per_g: list[tuple[float, Offer]] = []
    priced_raw: list[tuple[float, Offer]] = []
    for off in offers:
        if off.price is None or off.currency != currency:
            continue
        ppg = offer_price_per_gram(off, currency=currency)
        if ppg is not None:
            priced_per_g.append((ppg, off))
        else:
            priced_raw.append((off.price, off))

    if priced_per_g:
        ppg, off = min(priced_per_g, key=lambda t: t[0])
        pack_str = off.pack_size.raw_unit_str if off.pack_size else ""
        return Quantity(
            value=round(ppg, 4),
            unit=f"{currency}/g",
            raw_unit_str=f"{off.price} {off.currency} / {pack_str}".strip(),
        )
    if priced_raw:
        # Raw pack prices are only comparable when every offer has the same known pack size;
        # a lone offer is returned as-is. Otherwise there is no meaningful "cheapest".
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


def cheapest_offer_for_target(offers: list[Offer], currency: str = "USD") -> Optional[Quantity]:
    """Return cheapest listed price among offers, normalized to currency/g when mass pack_size is known.

    Falls back to raw pack price only when all such offers share one pack size. If no offer is
    comparable in `currency`, offers are compared only when they all use a single other currency;
    prices in mixed currencies are never compared and yield None. Returns None if all prices are null.
    """
    found = _cheapest_in_currency(offers, currency)
    if found is not None:
        return found
    other_currencies = {
        off.currency for off in offers if off.price is not None and off.currency is not None
    }
    if len(other_currencies) == 1:
        (only,) = other_currencies
        if only != currency:
            return _cheapest_in_currency(offers, only)
    return None


def estimate_route_cost(
    target_smiles: str,
    leaves: list[str],
    leaf_offers: dict[str, list[Offer]],
    currency: str = "USD",
) -> tuple[Optional[Quantity], bool]:
    """Estimate rough route cost per 1 g of target from real listed leaf prices only.

    Returns (est_cost, est_cost_complete). Never invents missing prices:
    - If 0 leaves have a known USD/g price -> (None, False)
    - If some leaves have prices and others do not -> (partial Quantity, False)
    - If all leaves have prices -> (Quantity, True)
    """
    est, complete, _ = estimate_route_cost_detailed(target_smiles, leaves, leaf_offers, currency)
    return est, complete


def estimate_route_cost_detailed(
    target_smiles: str,
    leaves: list[str],
    leaf_offers: dict[str, list[Offer]],
    currency: str = "USD",
) -> tuple[Optional[Quantity], bool, list[str]]:
    """Like estimate_route_cost, but also return one note per leaf left out of the sum and why.

    Notes have the form "<leaf SMILES> (unparseable SMILES)" or
    "<leaf SMILES> (no listed <currency>/g price)".
    """
    if not leaves:
        return None, False, []
    target_mol = normalize_smiles(target_smiles)
    if target_mol is None or target_mol.mol_weight <= 0:
        return None, False, [f"{target_smiles} (target unparseable; cost not estimated)"]

    counts = Counter(leaves)
    total = 0.0
    priced_leaves = 0
    total_unique_leaves = len(counts)
    skipped: list[str] = []

    for leaf_smi, count in counts.items():
        leaf_mol = normalize_smiles(leaf_smi)
        if leaf_mol is None:
            skipped.append(f"{leaf_smi} (unparseable SMILES)")
            continue
        offers = leaf_offers.get(leaf_mol.inchikey, [])
        ppg_candidates = [
            ppg for off in offers if (ppg := offer_price_per_gram(off, currency=currency)) is not None
        ]
        if not ppg_candidates:
            skipped.append(f"{leaf_smi} (no listed {currency}/g price)")
            continue
        best_ppg = min(ppg_candidates)
        equiv_mass_g = count * (leaf_mol.mol_weight / target_mol.mol_weight)
        total += best_ppg * equiv_mass_g
        priced_leaves += 1

    if priced_leaves == 0:
        return None, False, skipped

    complete = priced_leaves == total_unique_leaves
    suffix = "per 1 g target" if complete else "partial lower bound; per 1 g target"
    return (
        Quantity(
            value=round(total, 4),
            unit=currency,
            raw_unit_str=f"{round(total, 4)} {currency} ({suffix})",
        ),
        complete,
        skipped,
    )


SKIPPED_LEAVES_PREFIX = "; skipped leaves: "


def route_cost_method(skipped: list[str]) -> str:
    """Return EST_COST_METHOD, annotated with the leaves excluded from the estimate."""
    if not skipped:
        return EST_COST_METHOD
    return EST_COST_METHOD + SKIPPED_LEAVES_PREFIX + "; ".join(skipped)


def skipped_leaves_from_method(method: str) -> Optional[str]:
    """Inverse of route_cost_method: return the skipped-leaves note, or None if nothing was skipped."""
    if SKIPPED_LEAVES_PREFIX not in method:
        return None
    return method.split(SKIPPED_LEAVES_PREFIX, 1)[1]
