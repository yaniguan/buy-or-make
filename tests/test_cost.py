"""Tests for honest cost estimation: null prices stay null; partial leaves mark incomplete."""

from __future__ import annotations

from buyormake.canonical import normalize_smiles
from buyormake.cost import cheapest_offer_for_target, estimate_route_cost
from buyormake.schema import Offer, Quantity


def _offer(smiles: str, price=None, currency=None, pack_val=None, pack_unit="g") -> Offer:
    nm = normalize_smiles(smiles)
    assert nm is not None
    pack = (
        Quantity(value=pack_val, unit=pack_unit, raw_unit_str=f"{pack_val} {pack_unit}")
        if pack_val is not None
        else None
    )
    return Offer(
        inchikey=nm.inchikey,
        vendor="TestVendor",
        pack_size=pack,
        price=price,
        currency=currency,
        price_basis="listed_catalog" if price is not None else "unknown",
        source_id="local_stock",
        source_url="https://example.org/item",
        fetched_utc="2026-11-18T00:00:00Z",
        license="CC0-1.0",
    )


def test_price_null_stays_null() -> None:
    off = _offer("CCO", price=None, currency=None)
    assert off.price is None
    assert cheapest_offer_for_target([off]) is None


def test_cost_estimate_incomplete_when_leaf_has_no_price() -> None:
    target = "CC(=O)Oc1ccccc1C(=O)O"
    leaf1 = "CC(=O)O"
    leaf2 = "O=C(O)c1ccccc1O"
    nm1 = normalize_smiles(leaf1)
    nm2 = normalize_smiles(leaf2)
    assert nm1 and nm2

    # Case 1: neither leaf has a price
    est0, complete0 = estimate_route_cost(
        target,
        [leaf1, leaf2],
        {
            nm1.inchikey: [_offer(leaf1, price=None)],
            nm2.inchikey: [_offer(leaf2, price=None)],
        },
    )
    assert est0 is None
    assert complete0 is False

    # Case 2: leaf1 has price ($10 for 100 g = $0.10/g), leaf2 has price=None
    est1, complete1 = estimate_route_cost(
        target,
        [leaf1, leaf2],
        {
            nm1.inchikey: [_offer(leaf1, price=10.0, currency="USD", pack_val=100.0, pack_unit="g")],
            nm2.inchikey: [_offer(leaf2, price=None)],
        },
    )
    assert est1 is not None
    assert est1.value is not None and est1.value > 0
    assert complete1 is False

    # Case 3: both leaves have prices (leaf2 $5 for 500 mg = $10/g)
    est2, complete2 = estimate_route_cost(
        target,
        [leaf1, leaf2],
        {
            nm1.inchikey: [_offer(leaf1, price=10.0, currency="USD", pack_val=100.0, pack_unit="g")],
            nm2.inchikey: [_offer(leaf2, price=5.0, currency="USD", pack_val=500.0, pack_unit="mg")],
        },
    )
    assert est2 is not None
    assert complete2 is True
    assert est2.value > est1.value
