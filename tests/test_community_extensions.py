"""Unit and benchmark tests for buy-or-make community_extensions (Ideas #05 and #06)."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

_REPO_ROOT = str(Path(__file__).resolve().parents[1])
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

_ce = importlib.import_module("community_extensions")
AtomCountVec = _ce.AtomCountVec
GreenProcurementOption = _ce.GreenProcurementOption
SynthonPack = _ce.SynthonPack
TargetAnalogueSpec = _ce.TargetAnalogueSpec
compute_pareto_frontier = _ce.compute_pareto_frontier
optimize_independent_per_target_cart = _ce.optimize_independent_per_target_cart
optimize_shared_synthon_cart = _ce.optimize_shared_synthon_cart
pareto_dominates = _ce.pareto_dominates
run_pareto_green_benchmark = _ce.run_pareto_green_benchmark
run_plate_pack_benchmark = _ce.run_plate_pack_benchmark
select_pareto_green_option = _ce.select_pareto_green_option
select_scalar_cost_only = _ce.select_scalar_cost_only
verify_option_guardrails = _ce.verify_option_guardrails


def _make_valid_option(
    option_id: str,
    usd_per_g: float | None,
    lead_days: int,
    pmi: float,
    *,
    is_direct_buy: bool = False,
    bond_is_acyclic: bool = True,
    cost_complete: bool = True,
    target_c: int = 10,
) -> GreenProcurementOption:
    return GreenProcurementOption(
        option_id=option_id,
        is_direct_buy=is_direct_buy,
        steps=0 if is_direct_buy else 1,
        max_steps=2,
        bond_is_acyclic=bond_is_acyclic,
        target_atoms=AtomCountVec(target_c, 1, 1, 0, 0),
        leaving_delta=AtomCountVec(0, 0, 0, 0, 0)
        if is_direct_buy
        else AtomCountVec(0, 0, 1, 0, 0),
        reactant_atoms=()
        if is_direct_buy
        else (AtomCountVec(6, 0, 2, 0, 0), AtomCountVec(4, 1, 0, 0, 0)),
        all_leaves_in_stock=True,
        cost_complete=cost_complete,
        is_usd_per_gram=usd_per_g is not None,
        usd_per_g=usd_per_g,
        lead_time_days=lead_days,
        pmi_kg_per_kg=pmi,
    )


def test_pareto_guardrails_enforce_null_price_and_stoichiometry() -> None:
    opt_ok = _make_valid_option("ok", 30.0, 4, 18.0)
    assert verify_option_guardrails(opt_ok, budget_usd_per_g=50.0) == (True, None)
    assert opt_ok.e_factor_kg_per_kg == pytest.approx(17.0)

    # Null price -> strict_homogeneous_budget abstention
    opt_null = _make_valid_option(
        "null_p", None, 2, 10.0, is_direct_buy=True, cost_complete=False
    )
    assert verify_option_guardrails(opt_null, budget_usd_per_g=50.0) == (
        False,
        "strict_homogeneous_budget",
    )

    # Endocyclic ring cleavage -> acyclic_atom_balanced_route rejection
    opt_ring = _make_valid_option("ring_bad", 12.0, 3, 15.0, bond_is_acyclic=False)
    assert verify_option_guardrails(opt_ring, budget_usd_per_g=50.0) == (
        False,
        "acyclic_atom_balanced_route",
    )

    # Heavy-atom conservation mismatch -> acyclic_atom_balanced_route rejection
    opt_imbal = _make_valid_option("atom_bad", 12.0, 3, 15.0, target_c=11)
    assert verify_option_guardrails(opt_imbal, budget_usd_per_g=50.0) == (
        False,
        "acyclic_atom_balanced_route",
    )


def test_pareto_dominance_and_frontier_quarantine() -> None:
    fast_green = _make_valid_option("fast_green", 32.0, 4, 16.0)
    dominated_slow = _make_valid_option("dominated_slow", 35.0, 9, 55.0)
    cheap_dirty = _make_valid_option("cheap_dirty", 18.0, 14, 90.0)
    unpriced = _make_valid_option(
        "unpriced", None, 1, 8.0, is_direct_buy=True, cost_complete=False
    )
    ring_trap = _make_valid_option("ring_trap", 9.0, 2, 10.0, bond_is_acyclic=False)

    assert pareto_dominates(fast_green, dominated_slow) is True
    assert pareto_dominates(fast_green, fast_green) is False
    assert pareto_dominates(fast_green, cheap_dirty) is False

    res = compute_pareto_frontier(
        [fast_green, dominated_slow, cheap_dirty, unpriced, ring_trap],
        budget_usd_per_g=60.0,
    )
    frontier_ids = [o.option_id for o in res.pareto_frontier]
    assert frontier_ids == ["cheap_dirty", "fast_green"]
    assert res.abstained_unpriced_ids == ("unpriced",)
    assert res.rejected_stoichiometry_ids == ("ring_trap",)

    chosen = select_pareto_green_option(
        [fast_green, dominated_slow, cheap_dirty, unpriced, ring_trap],
        budget_usd_per_g=60.0,
        max_lead_days=10,
    )
    assert chosen.selected_option is not None
    assert chosen.selected_option.option_id == "fast_green"
    assert (
        select_scalar_cost_only([fast_green, cheap_dirty], budget_usd_per_g=60.0)
        == cheap_dirty
    )


def test_pareto_green_benchmark_metrics() -> None:
    summary = run_pareto_green_benchmark(
        n_candidates=80, budget_usd_per_g=60.0, max_lead_days=10, seed=42
    )
    assert summary["n_candidates"] == 80
    scalar_m = summary["scalar_cost_only"]
    pareto_m = summary["pareto_green_router"]
    assert isinstance(scalar_m, dict) and isinstance(pareto_m, dict)

    assert pareto_m["mean_pmi_waste_reduction_pct"] > 65.0
    assert pareto_m["on_time_delivery_rate"] == 1.0
    assert scalar_m["on_time_delivery_rate"] < 0.20
    assert pareto_m["within_budget_rate"] == 1.0
    assert pareto_m["null_price_leaks"] == 0
    assert pareto_m["stoichiometry_leaks"] == 0


def test_shared_synthon_cart_respects_budget_and_excludes_null_prices() -> None:
    synthons = (
        SynthonPack("A1", "amine", 50.0),
        SynthonPack("A2", "amine", 50.0),
        SynthonPack("B1", "acid", 50.0),
        SynthonPack("B2", "acid", 50.0),
        SynthonPack("EXP_1", "bespoke", 95.0),
        SynthonPack("EXP_2", "bespoke", 95.0),
        SynthonPack("NULL_BB", "acid", None),
    )
    targets = (
        TargetAnalogueSpec("T_exp", ("EXP_1", "EXP_2")),
        TargetAnalogueSpec("T_11", ("A1", "B1")),
        TargetAnalogueSpec("T_12", ("A1", "B2")),
        TargetAnalogueSpec("T_21", ("A2", "B1")),
        TargetAnalogueSpec("T_22", ("A2", "B2")),
        TargetAnalogueSpec("T_null", ("A1", "NULL_BB")),
    )

    indep = optimize_independent_per_target_cart(
        synthons, targets, cart_budget_usd=200.0
    )
    shared = optimize_shared_synthon_cart(synthons, targets, cart_budget_usd=200.0)

    # Independent greedy buyer spends $190 on T_exp (unlocking 1 target), whereas shared cart
    # buys A1, A2, B1, B2 ($200) and unlocks all 4 combinatorial analogues (and never T_null)!
    assert len(indep.unlocked_target_ids) == 1
    assert len(shared.unlocked_target_ids) == 4
    assert "T_null" not in shared.unlocked_target_ids
    assert shared.total_spend_usd <= 200.0


def test_plate_pack_procurement_benchmark_metrics() -> None:
    summary = run_plate_pack_benchmark(cart_budget_usd=500.0, seed=42)
    assert summary["cart_budget_usd"] == 500.0
    indep_m = summary["independent_per_target_buy"]
    shared_m = summary["shared_synthon_cart"]
    assert isinstance(indep_m, dict) and isinstance(shared_m, dict)

    assert shared_m["budget_respected"] is True
    assert shared_m["total_spend_usd"] <= 500.0
    assert (
        shared_m["target_analogues_per_100_usd"]
        >= 3.0 * indep_m["target_analogues_per_100_usd"]
    )
    assert shared_m["synthon_reuse_factor"] > 2.5
