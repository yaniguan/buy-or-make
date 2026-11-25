"""Tests for the Layer L5 Auto-Lab wire contract (`autolab_wire.py`)."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

from buyormake.routes import TemplateFallbackBackend
from buyormake.schema import Route

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

autolab_wire = importlib.import_module("autolab_wire")


def test_buyable_aspirin_triage_and_over_budget_block() -> None:
    assert autolab_wire.REPO_ID == "buy-or-make"
    assert autolab_wire.LAYER == "L5"
    assert autolab_wire.REPO_DIR == _REPO_ROOT

    raw_packet = {"goal_text": "Verify aspirin procurement within $25/g"}
    out = autolab_wire.process_packet(raw_packet)

    # Original packet must not be mutated
    assert raw_packet == {"goal_text": "Verify aspirin procurement within $25/g"}

    assert out["schema_version"] == "autolab-wire/v1"
    assert out["halted"] is False
    assert out["halt_reason"] is None
    proc = out["procurement"]
    assert proc["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"
    assert proc["status"] == "buyable"
    assert proc["listed_price_usd_per_g"] == pytest.approx(0.42)
    assert proc["n_steps"] == 0
    assert proc["price_basis"] == "listed_catalog"
    assert proc["within_budget"] is True
    assert proc["lead_time_days"] == 2
    assert proc["lead_time_complete"] is True
    assert proc["release_time_minutes"] == 2880
    assert proc["moq_pack_g"] == pytest.approx(1.0)
    assert proc["selected_pack_size"] is not None
    assert proc["selected_pack_size"]["value"] == pytest.approx(1.0)
    assert proc["selected_pack_size"]["unit"] == "g"
    assert proc["min_order_outlay_usd"] == pytest.approx(0.42)
    assert proc["precursor_smiles"] == []
    assert proc["precursor_basis"] == []
    assert proc["reaction_steps"] == []

    assert len(out["trace"]) == 1
    assert out["trace"][0]["repo"] == "buy-or-make"
    assert out["trace"][0]["layer"] == "L5"
    assert out["trace"][0]["status"] == "PASS"
    assert "release_min=2880" in out["trace"][0]["summary"]
    assert len(out["bottom_reasoning_traces"]) == 1
    assert out["bottom_reasoning_traces"][0]["repo"] == "buy-or-make"
    assert out["bottom_reasoning_traces"][0]["invariant"] == "unpriced_never_passes_budget"

    # Priced aspirin ($0.42/g) under strict $0.10/g budget must halt with L5_PROCUREMENT_BUDGET_BLOCK
    tight = autolab_wire.process_packet(
        {
            "candidate_smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "budget_usd_per_g": 0.10,
            "require_priced_within_budget": True,
        }
    )
    assert tight["halted"] is True
    assert tight["halt_reason"] == "L5_PROCUREMENT_BUDGET_BLOCK"
    assert tight["procurement"]["listed_price_usd_per_g"] == pytest.approx(0.42)
    assert tight["procurement"]["within_budget"] is False
    assert tight["trace"][-1]["status"] == "BLOCK"


def test_unlisted_smiles_null_price_honesty_and_budget_block() -> None:
    exotic_smiles = "C12C3C4C1C5C2C3C45"  # cubane: unlisted and not makeable from reference stock

    # Permissive inspection mode: honest null price, status="unknown", within_budget=False, not halted
    honest = autolab_wire.process_packet(
        {
            "candidate_smiles": exotic_smiles,
            "budget_usd_per_g": 1000.0,
            "require_priced_within_budget": False,
        }
    )
    assert honest["halted"] is False
    assert honest["halt_reason"] is None
    h_proc = honest["procurement"]
    assert h_proc["smiles"] == exotic_smiles
    assert h_proc["status"] == "unknown"
    assert h_proc["listed_price_usd_per_g"] is None
    assert h_proc["n_steps"] is None
    assert h_proc["price_basis"] == "unknown"
    assert h_proc["within_budget"] is False
    assert h_proc["lead_time_days"] is None
    assert h_proc["lead_time_complete"] is False
    assert h_proc["release_time_minutes"] is None
    assert h_proc["moq_pack_g"] is None
    assert h_proc["selected_pack_size"] is None
    assert h_proc["min_order_outlay_usd"] is None
    assert h_proc["precursor_smiles"] == []
    assert h_proc["precursor_basis"] == []
    assert h_proc["reaction_steps"] == []
    assert honest["trace"][-1]["status"] == "PASS"

    # Strict Lean budget mode (`unpriced_never_passes_budget`): null price blocks execution
    blocked = autolab_wire.process_packet(
        {
            "candidate_smiles": exotic_smiles,
            "budget_usd_per_g": 1000.0,
            "require_priced_within_budget": True,
        }
    )
    assert blocked["halted"] is True
    assert blocked["halt_reason"] == "L5_PROCUREMENT_BUDGET_BLOCK"
    assert blocked["procurement"]["status"] == "unknown"
    assert blocked["procurement"]["listed_price_usd_per_g"] is None
    assert blocked["procurement"]["price_basis"] == "unknown"
    assert blocked["procurement"]["within_budget"] is False
    assert blocked["trace"][-1]["status"] == "BLOCK"
    assert blocked["bottom_reasoning_traces"][-1]["invariant"] == "unpriced_never_passes_budget"


def test_makeable_acyclic_stoichiometric_route_from_reference_building_blocks() -> None:
    # Benzanilide: 1-step acyclic amide coupling from benzoic acid ($0.15/g, 1d) + aniline ($0.20/g, 1d)
    out = autolab_wire.process_packet(
        {
            "candidate_smiles": "O=C(Nc1ccccc1)c1ccccc1",
            "budget_usd_per_g": 5.0,
            "require_priced_within_budget": True,
        }
    )
    assert out["halted"] is False
    proc = out["procurement"]
    assert proc["status"] == "makeable"
    assert proc["n_steps"] == 1
    assert proc["price_basis"] == "listed_catalog"
    assert proc["listed_price_usd_per_g"] is not None
    assert 0.10 < proc["listed_price_usd_per_g"] < 0.50
    assert proc["within_budget"] is True
    assert proc["lead_time_days"] == 1
    assert proc["lead_time_complete"] is True
    assert proc["release_time_minutes"] == 1680  # 1 * 1440 + 1 * 240
    assert proc["moq_pack_g"] == pytest.approx(1.0)
    assert proc["min_order_outlay_usd"] == pytest.approx(0.35)
    assert set(proc["precursor_smiles"]) == {"O=C(O)c1ccccc1", "Nc1ccccc1"}
    assert len(proc["precursor_basis"]) == 2
    assert {p["catalog_id"] for p in proc["precursor_basis"]} == {"REF-BZA-007", "REF-ANI-008"}
    assert len(proc["reaction_steps"]) == 1
    assert proc["reaction_steps"][0]["leaving_delta"] == {"O": 1}


def test_two_step_suzuki_amide_propagates_leaf_lead_time_and_precursor_basis() -> None:
    # Issue #30 regression: 2-step Suzuki + amide target with 3-day bottleneck leaf (REF-PYB-011)
    pkt = autolab_wire.process_packet(
        {"candidate_smiles": "COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1", "max_steps": 2}
    )
    assert pkt["halted"] is False
    proc = pkt["procurement"]
    assert proc["status"] == "makeable"
    assert proc["n_steps"] == 2
    assert proc["lead_time_days"] == 3
    assert proc["lead_time_complete"] is True
    assert proc["release_time_minutes"] == 4800  # 3 * 1440 + 2 * 240
    assert "release_min=4800" in pkt["trace"][0]["summary"]
    assert proc["moq_pack_g"] == pytest.approx(1.0)
    assert proc["min_order_outlay_usd"] == pytest.approx(4.70)
    assert set(proc["precursor_smiles"]) == {
        "COc1ccc(C(=O)O)cc1",
        "Nc1ccc(Br)cc1",
        "OB(O)c1ccccn1",
    }
    assert {p["catalog_id"] for p in proc["precursor_basis"]} == {
        "REF-MBA-009",
        "REF-BRA-010",
        "REF-PYB-011",
    }
    deltas = [step["leaving_delta"] for step in proc["reaction_steps"]]
    assert {"O": 1} in deltas
    assert {"Br": 1, "B": 1, "O": 2} in deltas


def test_unimolecular_lactone_ring_opening_accepted_by_wire_validator() -> None:
    # Issue #27 regression: valid 1-step gamma-butyrolactone hydrolytic ring opening
    tb = TemplateFallbackBackend()
    gbl_steps = tb.one_step("O=C1CCCO1")
    assert len(gbl_steps) >= 1
    gbl_step = gbl_steps[0]
    assert gbl_step.reactants == ["O=C(O)CCCO"]
    route = Route(
        target_inchikey="YEJRWHAVMIAJKC-UHFFFAOYSA-N",
        n_steps=1,
        steps=[gbl_step],
        leaves=gbl_step.reactants,
        leaves_purchasable=[True],
        method="template_fallback",
        est_cost_method="leaf sum",
        est_cost_complete=True,
    )
    assert autolab_wire._route_is_acyclic_and_atom_balanced(route, max_steps=2) is True


def test_halted_packet_passthrough() -> None:
    upstream_halted = {
        "schema_version": "autolab-wire/v1",
        "candidate_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "halted": True,
        "halt_reason": "L6_ABSTAIN",
        "trace": [{"repo": "lab-decide", "layer": "L6", "status": "ABSTAIN", "summary": "OOD"}],
        "bottom_reasoning_traces": [],
    }
    out = autolab_wire.process_packet(upstream_halted)
    assert out == upstream_halted
    assert out is not upstream_halted
    assert "procurement" not in out


def test_cwd_independence_leaves_zero_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    out = autolab_wire.process_packet({"require_priced_within_budget": True})
    assert out["halted"] is False
    assert out["procurement"]["listed_price_usd_per_g"] == pytest.approx(0.42)
    assert list(tmp_path.iterdir()) == []
