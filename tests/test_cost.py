"""Tests for honest cost estimation: null prices stay null; partial leaves mark incomplete."""

from __future__ import annotations

from buyormake.canonical import normalize_smiles
from buyormake.cost import (
    EST_COST_METHOD,
    cheapest_offer_for_target,
    estimate_route_cost,
    estimate_route_cost_detailed,
    route_cost_method,
    skipped_leaves_from_method,
)
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


def test_cheapest_offer_never_compares_mixed_currencies() -> None:
    # No USD offer; JPY and EUR prices must not be compared as raw numbers.
    jpy = _offer("CCO", price=500.0, currency="JPY", pack_val=1.0, pack_unit="g")
    eur = _offer("CCO", price=30.0, currency="EUR", pack_val=1.0, pack_unit="g")
    assert cheapest_offer_for_target([jpy, eur]) is None

    # A single non-USD currency is compared within that currency only.
    eur_big = _offer("CCO", price=50.0, currency="EUR", pack_val=10.0, pack_unit="g")
    best = cheapest_offer_for_target([eur, eur_big])
    assert best is not None
    assert best.unit == "EUR/g"
    assert best.value == 5.0


def test_cheapest_offer_raw_pack_requires_same_pack_size() -> None:
    # Non-mass pack units cannot be normalized to USD/g.
    small = _offer("CCO", price=5.0, currency="USD", pack_val=1.0, pack_unit="mL")
    large = _offer("CCO", price=40.0, currency="USD", pack_val=1000.0, pack_unit="mL")
    assert cheapest_offer_for_target([small, large]) is None

    same_a = _offer("CCO", price=7.0, currency="USD", pack_val=1.0, pack_unit="mL")
    best = cheapest_offer_for_target([small, same_a])
    assert best is not None
    assert best.unit == "USD"
    assert best.value == 5.0


def test_cost_estimate_reports_skipped_leaves() -> None:
    target = "CC(=O)Oc1ccccc1C(=O)O"
    leaf1 = "CC(=O)O"
    leaf2 = "O=C(O)c1ccccc1O"
    bad_leaf = "not_a_smiles"
    nm1 = normalize_smiles(leaf1)
    nm2 = normalize_smiles(leaf2)
    assert nm1 and nm2

    est, complete, skipped = estimate_route_cost_detailed(
        target,
        [leaf1, leaf2, bad_leaf],
        {
            nm1.inchikey: [_offer(leaf1, price=10.0, currency="USD", pack_val=100.0, pack_unit="g")],
            nm2.inchikey: [_offer(leaf2, price=None)],
        },
    )
    assert est is not None
    assert complete is False
    assert skipped == [
        "O=C(O)c1ccccc1O (no listed USD/g price)",
        "not_a_smiles (unparseable SMILES)",
    ]

    method = route_cost_method(skipped)
    assert method.startswith(EST_COST_METHOD)
    assert skipped_leaves_from_method(method) == (
        "O=C(O)c1ccccc1O (no listed USD/g price); not_a_smiles (unparseable SMILES)"
    )
    assert route_cost_method([]) == EST_COST_METHOD
    assert skipped_leaves_from_method(EST_COST_METHOD) is None


def test_bulk_drum_moq_outlay_and_max_pack_multiple_gate() -> None:
    import pytest

    from buyormake.cost import estimate_route_cost_bounds
    from buyormake.rank import filter_decisions_from_ranked
    from buyormake.schema import RankResult

    target = "O=C(Nc1ccccc1)c1ccccc1"
    acid = "O=C(O)c1ccccc1"
    amine = "Nc1ccccc1"
    acid_ik = normalize_smiles(acid).inchikey  # type: ignore[union-attr]
    amine_ik = normalize_smiles(amine).inchikey  # type: ignore[union-attr]

    leaf_offers = {
        acid_ik: [
            _offer(acid, price=500.0, currency="USD", pack_val=1.0, pack_unit="kg"),
            _offer(acid, price=45.0, currency="USD", pack_val=1.0, pack_unit="g"),
        ],
        amine_ik: [
            _offer(amine, price=1.0, currency="USD", pack_val=1.0, pack_unit="g"),
        ],
    }

    # Uncapped unit price sees $0.50/g drum lower bound, while discrete MOQ outlay is $45 + $1 = $46
    bounds = estimate_route_cost_bounds(target, [acid, amine], leaf_offers)
    assert bounds.est_cost is not None and bounds.est_cost.value < 1.0
    assert bounds.min_outlay_cost is not None
    assert bounds.min_outlay_cost.value == pytest.approx(46.0)

    # Capping max_pack_multiple=50 on a 1 g synthesis excludes the 1 kg drum and uses the $45/1 g vial
    capped = estimate_route_cost_bounds(
        target, [acid, amine], leaf_offers, target_scale_g=1.0, max_pack_multiple=50.0
    )
    assert capped.est_cost is not None and capped.est_cost.value > 25.0

    # Budget filter distinguishes idealized 100%-yield `leaf_lower_bound` and blocks on MOQ outlay
    rr = RankResult(
        candidate=target,
        status="makeable",
        n_steps=1,
        est_cost=bounds.est_cost,
        est_cost_yield_adjusted=bounds.yield_adjusted_cost,
        min_outlay_cost=bounds.min_outlay_cost,
        step_yield_assumed=1.0,
        cost_complete=True,
    )
    lower_only = filter_decisions_from_ranked([rr], budget=10.0, require_complete_cost=True)
    assert lower_only.decisions[0].price_status == "leaf_lower_bound"

    moq_enforced = filter_decisions_from_ranked(
        [rr], budget=10.0, require_complete_cost=True, enforce_moq_outlay=True
    )
    assert moq_enforced.decisions[0].keep is False
    assert moq_enforced.decisions[0].over_budget is True
    assert any("over_budget_moq_outlay" in r for r in moq_enforced.decisions[0].reasons)


def test_step_depth_yield_compounding_penalizes_linear_over_convergent_routes() -> None:
    from buyormake.cost import estimate_route_cost_bounds
    from buyormake.rank import filter_decisions_from_ranked
    from buyormake.schema import RankResult, RouteStep

    target = "COc1ccc(C(=O)Nc2ccc(-c3ccc(NC(=O)C)cc3)cc2)cc1"
    l1 = "COc1ccc(C(=O)O)cc1"
    l2 = "Nc1ccc(Br)cc1"
    l3 = "OB(O)c1ccc(N)cc1"
    l4 = "CC(=O)O"
    leaves = [l1, l2, l3, l4]
    leaf_offers = {
        normalize_smiles(s).inchikey: [  # type: ignore[union-attr]
            _offer(s, price=10.0, currency="USD", pack_val=1.0, pack_unit="g")
        ]
        for s in leaves
    }

    # Linear 3-step route: l3 and l4 sit at depth d=3 (1 / 0.70^3 = 2.9155x mass!)
    linear_steps = [
        RouteStep(product=target, reactants=[l1, "Nc1ccc(-c2ccc(NC(=O)C)cc2)cc1"], template_or_reference="s1", source="t"),
        RouteStep(product="Nc1ccc(-c2ccc(NC(=O)C)cc2)cc1", reactants=[l2, "CC(=O)Nc1ccc(B(O)O)cc1"], template_or_reference="s2", source="t"),
        RouteStep(product="CC(=O)Nc1ccc(B(O)O)cc1", reactants=[l3, l4], template_or_reference="s3", source="t"),
    ]

    # Convergent 3-step route over the exact same 4 leaves: all 4 leaves sit at depth d=2 (1 / 0.70^2 = 2.0408x mass!)
    convergent_steps = [
        RouteStep(product=target, reactants=["COc1ccc(C(=O)Nc2ccc(Br)cc2)cc1", "CC(=O)Nc1ccc(B(O)O)cc1"], template_or_reference="s1", source="t"),
        RouteStep(product="COc1ccc(C(=O)Nc2ccc(Br)cc2)cc1", reactants=[l1, l2], template_or_reference="s2", source="t"),
        RouteStep(product="CC(=O)Nc1ccc(B(O)O)cc1", reactants=[l3, l4], template_or_reference="s3", source="t"),
    ]

    b_lin = estimate_route_cost_bounds(
        target, leaves, leaf_offers, steps=linear_steps, default_step_yield=0.70
    )
    b_conv = estimate_route_cost_bounds(
        target, leaves, leaf_offers, steps=convergent_steps, default_step_yield=0.70
    )
    assert b_lin.est_cost is not None and b_conv.est_cost is not None
    assert b_lin.lower_bound_cost == b_conv.lower_bound_cost  # Identical at 100% yield!
    assert b_lin.est_cost.value > b_conv.est_cost.value       # Linear strictly costlier at 70% step yield!
    assert b_lin.leaf_depths[l3] == [3]
    assert b_conv.leaf_depths[l3] == [2]

    rr_lin = RankResult(
        candidate=target,
        status="makeable",
        n_steps=3,
        est_cost=b_lin.est_cost,
        step_yield_assumed=0.70,
        cost_complete=True,
    )
    dec_ok = filter_decisions_from_ranked([rr_lin], budget=100.0, require_complete_cost=True)
    assert dec_ok.decisions[0].price_status == "verified"
