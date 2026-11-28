"""Tests for Lean 4.36.0 `BuyOrMakeLean` axiom audit and >=120-case Python differential parity."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SCRIPT_PATH = _REPO_ROOT / "scripts" / "check_lean_axioms_and_parity.py"


def _load_parity_module():
    spec = importlib.util.spec_from_file_location("check_lean_axioms_and_parity", _SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_lean_source_has_zero_sorry_and_zero_custom_axioms() -> None:
    mod = _load_parity_module()
    theorems = mod.extract_theorem_names()
    assert len(theorems) >= 18
    assert "checkBuyOrMake_ok_iff" in theorems
    assert "unpriced_never_passes_budget" in theorems
    assert "cross_halogen_swap_rejected" in theorems
    assert "suzuki_boron_imbalance_rejected" in theorems
    assert "shared_synthon_outlay_subadditive" in theorems
    assert "yield_compounded_monotone_depth" in theorems

    lean_readme = (_REPO_ROOT / "lean" / "README.md").read_text(encoding="utf-8")
    assert "no custom axioms" in lean_readme.lower()


def test_lean_lake_build_axioms_and_170_case_differential_parity() -> None:
    mod = _load_parity_module()
    if mod.find_lake_executable() is None:
        pytest.skip("Lean 4 `lake` toolchain not installed in this environment")

    summary = mod.run_lean_verification(verbose=False)
    assert summary["custom_axioms"] == 0
    assert summary["sorry_count"] == 0
    assert summary["mismatches"] == 0
    assert summary["differential_cases_total"] >= 170
    assert summary["suite_counts"]["shipping_gate"] >= 35
    assert summary["cases_against_shipping_code"] >= 80
    assert summary["cases_against_mirror"] >= 80
    assert (
        summary["cases_against_shipping_code"] + summary["cases_against_mirror"]
        == summary["differential_cases_total"]
    )
    assert set(summary["axioms_used_union"]).issubset(
        {"propext", "Classical.choice", "Quot.sound"}
    )
