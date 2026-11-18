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


def cheapest_offer_for_target(offers: list[Offer], currency: str = "USD") -> Optional[Quantity]:
    """Return cheapest listed price among offers (normalized to USD/g when mass pack_size is known,
    falling back to raw pack price if no mass-convertible offer exists). Returns None if all prices are null."""
    priced_per_g: list[tuple[float, Offer]] = []
    priced_raw: list[tuple[float, Offer]] = []
    for off in offers:
        if off.price is None or off.currency is None:
            continue
        ppg = offer_price_per_gram(off, currency=currency)
        if ppg is not None:
            priced_per_g.append((ppg, off))
        elif off.currency == currency:
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
        val, off = min(priced_raw, key=lambda t: t[0])
        pack_str = off.pack_size.raw_unit_str if (off.pack_size and off.pack_size.raw_unit_str) else "pack"
        return Quantity(
            value=round(val, 4),
            unit=currency,
            raw_unit_str=f"{val} {currency} / {pack_str}",
        )
    # Check if there is a non-USD priced offer
    any_priced = [off for off in offers if off.price is not None and off.currency is not None]
    if any_priced:
        off = min(any_priced, key=lambda o: o.price if o.price is not None else float("inf"))
        return Quantity(
            value=round(float(off.price), 4),
            unit=str(off.currency),
            raw_unit_str=f"{off.price} {off.currency}",
        )
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
    target_mol = normalize_smiles(target_smiles)
    if target_mol is None or target_mol.mol_weight <= 0 or not leaves:
        return None, False

    counts = Counter(leaves)
    total = 0.0
    priced_leaves = 0
    total_unique_leaves = len(counts)

    for leaf_smi, count in counts.items():
        leaf_mol = normalize_smiles(leaf_smi)
        if leaf_mol is None:
            continue
        offers = leaf_offers.get(leaf_mol.inchikey, [])
        ppg_candidates = [
            ppg for off in offers if (ppg := offer_price_per_gram(off, currency=currency)) is not None
        ]
        if not ppg_candidates:
            continue
        best_ppg = min(ppg_candidates)
        equiv_mass_g = count * (leaf_mol.mol_weight / target_mol.mol_weight)
        total += best_ppg * equiv_mass_g
        priced_leaves += 1

    if priced_leaves == 0:
        return None, False

    complete = priced_leaves == total_unique_leaves
    suffix = "per 1 g target" if complete else "partial lower bound; per 1 g target"
    return (
        Quantity(
            value=round(total, 4),
            unit=currency,
            raw_unit_str=f"{round(total, 4)} {currency} ({suffix})",
        ),
        complete,
    )
