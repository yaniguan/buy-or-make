#!/usr/bin/env python3
"""Verify `BuyOrMakeLean` (`Lean 4.36.0`) axioms and >=170-case differential parity against Python.

Checks:
1. Zero `sorry` or custom `axiom` declarations in `lean/BuyOrMakeLean.lean`.
2. `lake build` succeeds cleanly.
3. `#print axioms` on every theorem in `BuyOrMakeLean` uses ONLY standard Lean kernel axioms
   (subset of `{"propext", "Classical.choice", "Quot.sound"}`, zero `sorryAx`).
4. `>= 170` `#eval` differential parity test cases across:
   - Suite A (`proposals`, 60 cases — mirror + delegated `filter_decisions_from_ranked`):
     13-channel stoichiometric `AtomVec` conservation, acyclic bond rule, strict-null `USD/g` budget
     gating, and license policy (`checkBuyOrMake` vs `verify_option_guardrails` + RDKit),
   - Suite B (`moq`, 45 cases — shipping `buyormake.cost.min_outlay_for_offers`):
     Discrete pack-size MOQ outlay, multi-SKU selection, shared-synthon subadditivity, and MOQ budget
     gating (`minSkuOutlayCents` / `discretePackOutlayCents` / `moqOutlayPolicyOK`),
   - Suite C (`yield_pareto`, 30 cases — mirror `pareto_dominates` & fixed-point yield):
     Step-depth yield compounding & 3D Green Pareto dominance (`compoundYieldMassUg` & `paretoDominates`),
   - Suite D (`shipping_gate`, 40 cases — shipping `evaluate_candidates` + `filter_decisions_from_ranked`):
     End-to-end local-catalog triage and default budget gating over priced/unpriced/incomplete/non-USD-per-g/
     out-of-stock/no-budget cases compared directly against `BuyOrMakeLean.checkBuyOrMake`.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from buyormake.cost import min_outlay_for_offers  # noqa: E402
from buyormake.rank import evaluate_candidates, filter_decisions_from_ranked  # noqa: E402
from buyormake.routes import (  # noqa: E402
    RETRO_TEMPLATE_DELTAS,
    UNIMOLECULAR_RING_OPENING_TEMPLATES,
    TemplateFallbackBackend,
)
from buyormake.schema import Offer, Quantity  # noqa: E402
from community_extensions.pareto_green_triage import (  # noqa: E402
    AtomCountVec,
    GreenProcurementOption,
    atom_vec_from_counts,
    atom_vec_from_smiles,
    pareto_dominates,
    verify_option_guardrails,
)

LEAN_DIR = REPO_ROOT / "lean"
LEAN_FILE = LEAN_DIR / "BuyOrMakeLean.lean"
ALLOWED_AXIOMS: frozenset[str] = frozenset({"propext", "Classical.choice", "Quot.sound"})


def find_lake_executable() -> Optional[str]:
    """Locate `lake` on PATH or under `~/.elan/bin`."""
    found = shutil.which("lake")
    if found:
        return found
    home = Path.home()
    for cand in (home / ".elan" / "bin" / "lake.exe", home / ".elan" / "bin" / "lake"):
        if cand.exists():
            return str(cand)
    return None


def _strip_lean_comments(source: str) -> str:
    # Strip nested/block `/- ... -/` comments and `-- ...` line comments
    text = re.sub(r"/-.*?-/", " ", source, flags=re.DOTALL)
    text = re.sub(r"--[^\n]*", " ", text)
    return text


def extract_theorem_names(lean_path: Path = LEAN_FILE) -> list[str]:
    raw = lean_path.read_text(encoding="utf-8")
    code = _strip_lean_comments(raw)
    if re.search(r"\bsorry\b", code):
        raise AssertionError(f"Forbidden `sorry` token detected in {lean_path}")
    if re.search(r"(?m)^\s*axiom\s+", code):
        raise AssertionError(f"Forbidden custom `axiom` declaration detected in {lean_path}")
    names = re.findall(r"\btheorem\s+([A-Za-z0-9_']+)\b", code)
    if not names:
        raise AssertionError(f"No theorems found in {lean_path}")
    return list(dict.fromkeys(names))


def _atom_vec_lean_lit(v: AtomCountVec) -> str:
    return (
        f"{{ c := {v.c}, n := {v.n}, o := {v.o}, s := {v.s}, "
        f"f := {v.f}, cl := {v.cl}, br := {v.br}, i := {v.i}, "
        f"b := {v.b}, p := {v.p}, si := {v.si}, "
        f"metal := {v.metal}, otherHeavy := {v.other_heavy} }}"
    )


def _opt_nat_lean(val: Optional[int]) -> str:
    return "none" if val is None else f"(some {int(val)})"


def _bool_lean(val: bool) -> str:
    return "true" if val else "false"


def _python_check_proposal(case: dict[str, Any]) -> Optional[str]:
    """Exact Python oracle corresponding to `BuyOrMakeLean.checkBuyOrMake`."""
    budget_cents: Optional[int] = case["budget_cents"]
    price_cents: Optional[int] = case["price_cents"]
    is_direct_buy: bool = case["is_direct_buy"]
    steps: int = case["steps"]
    max_steps: int = case["max_steps"]
    bond_is_acyclic: bool = case["bond_is_acyclic"]
    all_leaves_in_stock: bool = case["all_leaves_in_stock"]
    cost_complete: bool = case["cost_complete"]
    is_usd_per_gram: bool = case["is_usd_per_gram"]
    target_atoms: AtomCountVec = case["target_atoms"]
    leaving_delta: AtomCountVec = case["leaving_delta"]
    reactant_atoms: tuple[AtomCountVec, ...] = case["reactant_atoms"]
    licenses: list[str] = case["licenses"]
    commercial_mode: bool = case["commercial_mode"]

    # 1. Stoichiometry check (cross-checked with community_extensions.verify_option_guardrails)
    dummy_opt = GreenProcurementOption(
        option_id=case["case_id"],
        is_direct_buy=is_direct_buy,
        steps=steps,
        max_steps=max_steps,
        bond_is_acyclic=bond_is_acyclic,
        target_atoms=target_atoms,
        leaving_delta=leaving_delta,
        reactant_atoms=reactant_atoms,
        all_leaves_in_stock=all_leaves_in_stock,
        cost_complete=True,
        is_usd_per_gram=True,
        usd_per_g=1.0,
        lead_time_days=1,
        pmi_kg_per_kg=5.0,
    )
    stoich_ok, stoich_rule = verify_option_guardrails(dummy_opt, budget_usd_per_g=None)
    if not stoich_ok:
        assert stoich_rule == "acyclic_atom_balanced_route"
        return "acyclic_atom_balanced_route"

    # 2. Strict-null homogeneous USD/g budget check
    if budget_cents is not None:
        if price_cents is None or not cost_complete or not is_usd_per_gram or price_cents > budget_cents:
            return "strict_homogeneous_budget"

    # 3. Transitive license gate
    if commercial_mode:
        lic_ok = all(lic == "permissive" for lic in licenses)
    else:
        lic_ok = all(lic != "evalOnly" for lic in licenses)
    if not lic_ok:
        return "transitive_license_gate"

    return None


def build_differential_cases(seed: int = 20261127) -> dict[str, list[dict[str, Any]]]:
    rng = random.Random(seed)
    tb = TemplateFallbackBackend()

    # ---------- Suite A: 60 Stoichiometric / Budget / License Proposal Cases ----------
    proposal_cases: list[dict[str, Any]] = []

    seed_molecules = [
        "O=C(Nc1ccccc1)c1ccccc1",
        "CC(=O)Oc1ccccc1C(=O)O",
        "Fc1ccc(-c2ccccc2)cc1",
        "COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1",
        "O=S(=O)(Nc1ccc(Cl)cc1)c1ccc(F)cc1",
        "CCOc1ccc(C(=O)OCC)cc1",
        "O=C(Nc1ccc(F)cc1)Nc1cccc(Br)c1",
        "CCOC(=O)Nc1ccc(OC)cc1",
        "c1ccc(CNc2ccc(Cl)cc2)cc1",
        "O=C1CCCO1",
        "O=C1CCCN1",
        "O=C1NC2(CC1)CCN(c1ncccn1)CC2",
    ]
    for smi in seed_molecules:
        t_vec = atom_vec_from_smiles(smi)
        for st in tb.one_step(smi)[:2]:
            tpl_name = st.template_or_reference.split(":", 1)[0].strip()
            d_vec = atom_vec_from_counts(RETRO_TEMPLATE_DELTAS[tpl_name])
            r_vecs = tuple(atom_vec_from_smiles(r) for r in st.reactants)
            proposal_cases.append(
                {
                    "case_id": f"rdkit_valid_{len(proposal_cases)}",
                    "is_direct_buy": False,
                    "steps": 1,
                    "max_steps": 3,
                    "bond_is_acyclic": True,
                    "target_atoms": t_vec,
                    "leaving_delta": d_vec,
                    "reactant_atoms": r_vecs,
                    "all_leaves_in_stock": True,
                    "cost_complete": True,
                    "is_usd_per_gram": True,
                    "price_cents": 850,
                    "budget_cents": 2500,
                    "licenses": ["permissive"],
                    "commercial_mode": True,
                }
            )

    # Adversarial halogen swaps, Boron omissions, Si/P/metal synthons, ring splits, budget & license traps
    suzuki_t = atom_vec_from_smiles("Fc1ccc(-c2ccccc2)cc1")
    suzuki_d = atom_vec_from_counts(RETRO_TEMPLATE_DELTAS["suzuki_biaryl"])
    adversarial_specs = [
        # Cross-halogen F -> Cl swap preserving total hal count
        {
            "is_direct_buy": False,
            "steps": 1,
            "max_steps": 2,
            "bond_is_acyclic": True,
            "target_atoms": suzuki_t,
            "leaving_delta": suzuki_d,
            "reactant_atoms": (atom_vec_from_smiles("Clc1ccc(Br)cc1"), atom_vec_from_smiles("OB(O)c1ccccc1")),
            "all_leaves_in_stock": True,
            "cost_complete": True,
            "is_usd_per_gram": True,
            "price_cents": 500,
            "budget_cents": 2500,
            "licenses": ["permissive"],
            "commercial_mode": True,
        },
        # Missing Boron in Suzuki coupling
        {
            "is_direct_buy": False,
            "steps": 1,
            "max_steps": 2,
            "bond_is_acyclic": True,
            "target_atoms": suzuki_t,
            "leaving_delta": suzuki_d,
            "reactant_atoms": (atom_vec_from_smiles("Fc1ccc(Br)cc1"), AtomCountVec(c=6, o=2, b=0)),
            "all_leaves_in_stock": True,
            "cost_complete": True,
            "is_usd_per_gram": True,
            "price_cents": 500,
            "budget_cents": 2500,
            "licenses": ["permissive"],
            "commercial_mode": True,
        },
        # Organosilicon / phosphonate / metal synthon mismatch
        {
            "is_direct_buy": False,
            "steps": 1,
            "max_steps": 2,
            "bond_is_acyclic": True,
            "target_atoms": atom_vec_from_smiles("C[Si](C)(C)Oc1ccccc1"),
            "leaving_delta": AtomCountVec(cl=1),
            "reactant_atoms": (AtomCountVec(c=3, cl=1, si=0), atom_vec_from_smiles("Oc1ccccc1")),
            "all_leaves_in_stock": True,
            "cost_complete": True,
            "is_usd_per_gram": True,
            "price_cents": 400,
            "budget_cents": 2500,
            "licenses": ["permissive"],
            "commercial_mode": True,
        },
    ]
    for spec in adversarial_specs:
        proposal_cases.append({"case_id": f"adv_{len(proposal_cases)}", **spec})

    license_pool = ["permissive", "shareAlike", "nonCommercial", "evalOnly", "customTerms"]
    while len(proposal_cases) < 60:
        is_buy = rng.random() < 0.25
        steps = 0 if is_buy else rng.randint(0, 6)
        max_steps = rng.randint(1, 5)
        t_vec = AtomCountVec(
            c=rng.randint(4, 18),
            n=rng.randint(0, 4),
            o=rng.randint(0, 4),
            s=rng.randint(0, 2),
            f=rng.randint(0, 3),
            cl=rng.randint(0, 2),
            br=rng.randint(0, 1),
            i=rng.randint(0, 1),
            b=rng.randint(0, 1),
            p=rng.randint(0, 1),
            si=rng.randint(0, 1),
            metal=rng.randint(0, 1),
        )
        d_vec = AtomCountVec(
            o=rng.randint(0, 2),
            cl=rng.randint(0, 1),
            br=rng.randint(0, 1),
            b=rng.randint(0, 1),
        )
        sum_vec = t_vec.add(d_vec)
        # Split sum_vec into two reactants, then optionally inject an atom/halogen/boron perturbation
        r1 = AtomCountVec(
            c=sum_vec.c // 2,
            n=sum_vec.n // 2,
            o=sum_vec.o // 2,
            s=sum_vec.s,
            f=sum_vec.f,
            cl=sum_vec.cl,
            b=sum_vec.b,
            si=sum_vec.si,
        )
        r2 = AtomCountVec(
            c=sum_vec.c - r1.c,
            n=sum_vec.n - r1.n,
            o=sum_vec.o - r1.o,
            br=sum_vec.br,
            i=sum_vec.i,
            p=sum_vec.p,
            metal=sum_vec.metal,
        )
        mode = rng.choice(["balanced", "swap_hal", "drop_b", "ring_bad", "null_price", "bad_unit", "over_bud", "lic"])
        if mode == "swap_hal":
            r1 = AtomCountVec(**{**r1.__dict__, "f": r1.f + 1})
            r2 = AtomCountVec(**{**r2.__dict__, "cl": max(0, r2.cl - 1)})
        elif mode == "drop_b":
            r1 = AtomCountVec(**{**r1.__dict__, "b": r1.b + 1})
        bond_acyclic = mode != "ring_bad"
        all_stock = rng.random() > 0.1
        price_cents: Optional[int] = None if mode == "null_price" else (4000 if mode == "over_bud" else 900)
        budget_cents: Optional[int] = None if rng.random() < 0.15 else 2500
        cost_complete = rng.random() > 0.15
        is_usdg = mode != "bad_unit"
        lics = [rng.choice(license_pool) for _ in range(rng.randint(1, 2))] if mode == "lic" else ["permissive"]
        comm_mode = rng.random() > 0.3

        proposal_cases.append(
            {
                "case_id": f" fuzz_prop_{len(proposal_cases)}",
                "is_direct_buy": is_buy,
                "steps": steps,
                "max_steps": max_steps,
                "bond_is_acyclic": bond_acyclic,
                "target_atoms": t_vec,
                "leaving_delta": AtomCountVec() if is_buy else d_vec,
                "reactant_atoms": () if is_buy else (r1, r2),
                "all_leaves_in_stock": all_stock,
                "cost_complete": cost_complete,
                "is_usd_per_gram": is_usdg,
                "price_cents": price_cents,
                "budget_cents": budget_cents,
                "licenses": lics,
                "commercial_mode": comm_mode,
            }
        )

    # ---------- Suite B: 45 Discrete MOQ Pack Outlay & Shared-Synthon Cases ----------
    moq_cases: list[dict[str, Any]] = [
        # Canonical Issue #28 witness: 0.62 g benzoic acid ($500/1kg drum vs $45/1g vial)
        {
            "case_id": "moq_drum_vs_vial",
            "req_ug": 620_000,
            "req_b_ug": 310_000,
            "skus": [(1_000_000_000, 50_000), (1_000_000, 4_500)],
            "max_outlay_cents": 1_000,
            "require_complete": True,
        },
        # Exact pack boundary: 1.000000 g vs 1.000001 g
        {
            "case_id": "moq_exact_1g",
            "req_ug": 1_000_000,
            "req_b_ug": 500_000,
            "skus": [(1_000_000, 1_200), (5_000_000, 4_000)],
            "max_outlay_cents": 2_000,
            "require_complete": True,
        },
        {
            "case_id": "moq_just_over_1g",
            "req_ug": 1_000_005,
            "req_b_ug": 500_000,
            "skus": [(1_000_000, 1_200), (5_000_000, 4_000)],
            "max_outlay_cents": 2_500,
            "require_complete": True,
        },
    ]
    pack_choices_ug = [100_000, 250_000, 500_000, 1_000_000, 5_000_000, 25_000_000, 1_000_000_000]
    while len(moq_cases) < 45:
        n_skus = rng.randint(0, 4)
        chosen_packs = rng.sample(pack_choices_ug, k=n_skus)
        skus = [(p, rng.randint(150, 60_000)) for p in chosen_packs]
        req_ug = rng.randint(25_000, 8_500_000)
        req_b_ug = rng.randint(25_000, 3_500_000)
        max_outlay = None if rng.random() < 0.2 else rng.randint(500, 25_000)
        moq_cases.append(
            {
                "case_id": f"moq_{len(moq_cases)}",
                "req_ug": req_ug,
                "req_b_ug": req_b_ug,
                "skus": skus,
                "max_outlay_cents": max_outlay,
                "require_complete": rng.choice([True, False]),
            }
        )

    # ---------- Suite C: 30 Yield-Compounding & 3D Pareto Dominance Cases ----------
    yield_pareto_cases: list[dict[str, Any]] = []
    for idx in range(30):
        base_ug = rng.randint(100_000, 2_000_000)
        yield_pm = rng.choice([500, 650, 700, 750, 800, 850, 900, 1000])
        depth = rng.randint(0, 5)
        a_cost = rng.randint(100, 5000)
        b_cost = a_cost if rng.random() < 0.3 else rng.randint(100, 5000)
        a_lead = rng.randint(1, 14)
        b_lead = a_lead if rng.random() < 0.3 else rng.randint(1, 14)
        a_pmi = rng.randint(5000, 90000)
        b_pmi = a_pmi if rng.random() < 0.3 else rng.randint(5000, 90000)
        yield_pareto_cases.append(
            {
                "case_id": f"yp_{idx}",
                "base_ug": base_ug,
                "yield_pm": yield_pm,
                "depth": depth,
                "a": (a_cost, a_lead, a_pmi),
                "b": (b_cost, b_lead, b_pmi),
            }
        )

    return {
        "proposals": proposal_cases,
        "moq": moq_cases,
        "yield_pareto": yield_pareto_cases,
        "shipping_gate": _build_shipping_gate_cases(),
    }


def _row(
    smiles: str,
    *,
    vendor: str = "LocalCat",
    pack_val: Optional[float] = 1.0,
    pack_unit: Optional[str] = "g",
    price: Optional[float] = 5.0,
    currency: Optional[str] = "USD",
    in_stock: Optional[bool] = True,
) -> dict[str, Any]:
    return {
        "smiles": smiles,
        "vendor": vendor,
        "pack_val": pack_val,
        "pack_unit": pack_unit if pack_val is not None else None,
        "price": price,
        "currency": currency if price is not None else None,
        "in_stock": in_stock,
    }


def _build_shipping_gate_cases() -> list[dict[str, Any]]:
    """Build 40 end-to-end shipping-code cases (`evaluate_candidates` + `filter_decisions_from_ranked`)."""
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    paracetamol = "CC(=O)Nc1ccc(O)cc1"
    ibuprofen = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
    caffeine = "Cn1c(=O)c2c(ncn2C)n(C)c1=O"
    benzanilide = "O=C(Nc1ccccc1)c1ccccc1"
    benzoic_acid = "O=C(O)c1ccccc1"
    aniline = "Nc1ccccc1"
    fluorobiphenyl = "Fc1ccc(-c2ccccc2)cc1"
    bromofluorobenzene = "Fc1ccc(Br)cc1"
    phenylboronic_acid = "OB(O)c1ccccc1"
    phenetole = "CCOc1ccccc1"
    phenol = "Oc1ccccc1"
    bromoethane = "CCBr"
    sulfonamide = "O=S(=O)(Nc1ccc(Cl)cc1)c1ccc(F)cc1"
    f_sulfonyl_cl = "O=S(=O)(Cl)c1ccc(F)cc1"
    cl_aniline = "Nc1ccc(Cl)cc1"
    biaryl_amide = "COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1"
    mba = "COc1ccc(C(=O)O)cc1"
    bromoaniline = "Nc1ccc(Br)cc1"
    py_boronic = "OB(O)c1ccccn1"
    lactone = "O=C1CCCO1"
    hydroxybutyric = "O=C(O)CCCO"
    cubane = "C12C3C4C1C5C2C3C45"

    cases: list[dict[str, Any]] = [
        # --- 1..8: Direct buyable priced in USD/g under and at budget ---
        {"case_id": "ship_buy_under_1", "smiles": aspirin, "max_steps": 2, "budget": 25.0, "rows": [_row(aspirin, price=4.50)]},
        {"case_id": "ship_buy_under_mg_pack", "smiles": paracetamol, "max_steps": 2, "budget": 25.0, "rows": [_row(paracetamol, pack_val=500.0, pack_unit="mg", price=6.00)]},
        {"case_id": "ship_buy_under_kg_pack", "smiles": ibuprofen, "max_steps": 2, "budget": 5.0, "rows": [_row(ibuprofen, pack_val=1.0, pack_unit="kg", price=1200.0)]},
        {"case_id": "ship_buy_exact_boundary", "smiles": caffeine, "max_steps": 2, "budget": 15.0, "rows": [_row(caffeine, price=15.00)]},
        {"case_id": "ship_buy_in_stock_none_under", "smiles": aspirin, "max_steps": 1, "budget": 10.0, "rows": [_row(aspirin, price=8.25, in_stock=None)]},
        {"case_id": "ship_buy_mixed_oos_and_instock", "smiles": paracetamol, "max_steps": 1, "budget": 20.0, "rows": [_row(paracetamol, vendor="OOSV", price=1.00, in_stock=False), _row(paracetamol, vendor="OKV", price=14.00, in_stock=True)]},
        {"case_id": "ship_buy_multi_vendor_cheapest", "smiles": caffeine, "max_steps": 1, "budget": 12.0, "rows": [_row(caffeine, vendor="ExpV", price=30.00), _row(caffeine, vendor="CheapV", price=9.50)]},
        {"case_id": "ship_buy_salt_equiv_under", "smiles": "CCN(CC)CC(=O)Nc1c(C)cccc1C.Cl", "max_steps": 1, "budget": 25.0, "rows": [_row("CCN(CC)CC(=O)Nc1c(C)cccc1C", price=10.00)]},
        # --- 9..13: Direct buyable priced in USD/g OVER budget ---
        {"case_id": "ship_buy_over_1", "smiles": aspirin, "max_steps": 2, "budget": 10.0, "rows": [_row(aspirin, price=25.00)]},
        {"case_id": "ship_buy_over_mg_pack", "smiles": paracetamol, "max_steps": 2, "budget": 20.0, "rows": [_row(paracetamol, pack_val=100.0, pack_unit="mg", price=5.00)]},
        {"case_id": "ship_buy_over_just_above", "smiles": ibuprofen, "max_steps": 1, "budget": 10.0, "rows": [_row(ibuprofen, price=10.50)]},
        {"case_id": "ship_buy_mixed_oos_cheap_instock_over", "smiles": caffeine, "max_steps": 1, "budget": 15.0, "rows": [_row(caffeine, vendor="OOSCheap", price=2.00, in_stock=False), _row(caffeine, vendor="InStockExp", price=40.00, in_stock=True)]},
        {"case_id": "ship_buy_in_stock_none_over", "smiles": aspirin, "max_steps": 1, "budget": 1.0, "rows": [_row(aspirin, price=4.00, in_stock=None)]},
        # --- 14..18: Unpriced (`price=None`) under active budget (tests default require_complete_cost=True!) ---
        {"case_id": "ship_buy_null_price_aspirin", "smiles": aspirin, "max_steps": 2, "budget": 25.0, "rows": [_row(aspirin, pack_val=None, price=None)]},
        {"case_id": "ship_buy_null_price_with_pack", "smiles": paracetamol, "max_steps": 2, "budget": 100.0, "rows": [_row(paracetamol, pack_val=5.0, pack_unit="g", price=None)]},
        {"case_id": "ship_buy_null_price_in_stock_none", "smiles": caffeine, "max_steps": 1, "budget": 50.0, "rows": [_row(caffeine, pack_val=None, price=None, in_stock=None)]},
        {"case_id": "ship_make_all_leaves_unpriced", "smiles": benzanilide, "max_steps": 2, "budget": 25.0, "rows": [_row(benzoic_acid, pack_val=None, price=None), _row(aniline, pack_val=None, price=None)]},
        {"case_id": "ship_make_suzuki_all_leaves_unpriced", "smiles": fluorobiphenyl, "max_steps": 2, "budget": 50.0, "rows": [_row(bromofluorobenzene, pack_val=None, price=None), _row(phenylboronic_acid, pack_val=None, price=None)]},
        # --- 19..25: Complete-cost makeable routes (1-step & 2-step) under and over budget ---
        {"case_id": "ship_make_1step_amide_under", "smiles": benzanilide, "max_steps": 2, "budget": 10.0, "rows": [_row(benzoic_acid, price=1.00), _row(aniline, price=2.00)]},
        {"case_id": "ship_make_1step_amide_over", "smiles": benzanilide, "max_steps": 2, "budget": 5.0, "rows": [_row(benzoic_acid, price=10.00), _row(aniline, price=10.00)]},
        {"case_id": "ship_make_1step_suzuki_under", "smiles": fluorobiphenyl, "max_steps": 2, "budget": 20.0, "rows": [_row(bromofluorobenzene, price=3.00), _row(phenylboronic_acid, price=4.00)]},
        {"case_id": "ship_make_1step_williamson_under", "smiles": phenetole, "max_steps": 2, "budget": 15.0, "rows": [_row(phenol, price=2.00), _row(bromoethane, price=3.00)]},
        {"case_id": "ship_make_1step_sulfonamide_over", "smiles": sulfonamide, "max_steps": 2, "budget": 4.0, "rows": [_row(f_sulfonyl_cl, price=8.00), _row(cl_aniline, price=6.00)]},
        {"case_id": "ship_make_1step_lactone_under", "smiles": lactone, "max_steps": 2, "budget": 12.0, "rows": [_row(hydroxybutyric, price=5.00)]},
        {"case_id": "ship_make_2step_biaryl_amide_under", "smiles": biaryl_amide, "max_steps": 2, "budget": 25.0, "rows": [_row(mba, price=2.00), _row(bromoaniline, price=3.00), _row(py_boronic, price=6.00)]},
        {"case_id": "ship_make_2step_biaryl_amide_over", "smiles": biaryl_amide, "max_steps": 2, "budget": 5.0, "rows": [_row(mba, price=5.00), _row(bromoaniline, price=8.00), _row(py_boronic, price=12.00)]},
        # --- 27..30: Incomplete route cost (`cost_complete=False`) under active budget ---
        {"case_id": "ship_make_incomplete_amide", "smiles": benzanilide, "max_steps": 2, "budget": 25.0, "rows": [_row(benzoic_acid, price=1.00), _row(aniline, pack_val=None, price=None)]},
        {"case_id": "ship_make_incomplete_suzuki", "smiles": fluorobiphenyl, "max_steps": 2, "budget": 25.0, "rows": [_row(bromofluorobenzene, price=2.00), _row(phenylboronic_acid, pack_val=None, price=None)]},
        {"case_id": "ship_make_incomplete_2step", "smiles": biaryl_amide, "max_steps": 2, "budget": 50.0, "rows": [_row(mba, price=1.00), _row(bromoaniline, price=2.00), _row(py_boronic, pack_val=None, price=None)]},
        {"case_id": "ship_make_incomplete_raw_ml_leaf", "smiles": phenetole, "max_steps": 2, "budget": 25.0, "rows": [_row(phenol, price=2.00), _row(bromoethane, pack_val=10.0, pack_unit="mL", price=4.00)]},
        # --- 31..34: Non-`USD/g` unit under active USD budget ---
        {"case_id": "ship_buy_raw_ml_pack_usd", "smiles": "CCO", "max_steps": 1, "budget": 25.0, "rows": [_row("CCO", pack_val=1.0, pack_unit="mL", price=5.00, currency="USD")]},
        {"case_id": "ship_buy_eur_per_g", "smiles": aspirin, "max_steps": 1, "budget": 25.0, "rows": [_row(aspirin, pack_val=1.0, pack_unit="g", price=6.00, currency="EUR")]},
        {"case_id": "ship_buy_jpy_per_g", "smiles": caffeine, "max_steps": 1, "budget": 25.0, "rows": [_row(caffeine, pack_val=1.0, pack_unit="g", price=500.0, currency="JPY")]},
        {"case_id": "ship_buy_raw_bottle_pack", "smiles": ibuprofen, "max_steps": 1, "budget": 50.0, "rows": [_row(ibuprofen, pack_val=1.0, pack_unit="bottle", price=12.00, currency="USD")]},
        # --- 35..37: Out-of-stock (`in_stock=False`) & unreachable molecules ---
        {"case_id": "ship_oos_direct_offer", "smiles": aspirin, "max_steps": 2, "budget": 1000.0, "rows": [_row(aspirin, price=10.00, in_stock=False)]},
        {"case_id": "ship_oos_route_leaves", "smiles": benzanilide, "max_steps": 2, "budget": 100.0, "rows": [_row(benzoic_acid, price=1.00, in_stock=False), _row(aniline, price=1.00, in_stock=False)]},
        {"case_id": "ship_unlisted_cubane", "smiles": cubane, "max_steps": 2, "budget": 100.0, "rows": [_row(aspirin, price=1.00)]},
        # --- 38..40: No active budget (`budget=None`) ---
        {"case_id": "ship_nobudget_unpriced_buy", "smiles": aspirin, "max_steps": 1, "budget": None, "rows": [_row(aspirin, pack_val=None, price=None)]},
        {"case_id": "ship_nobudget_incomplete_make", "smiles": benzanilide, "max_steps": 2, "budget": None, "rows": [_row(benzoic_acid, price=1.00), _row(aniline, pack_val=None, price=None)]},
        {"case_id": "ship_nobudget_oos_still_rejected", "smiles": paracetamol, "max_steps": 1, "budget": None, "rows": [_row(paracetamol, price=5.00, in_stock=False)]},
    ]
    return cases


def _eval_shipping_gate_case(
    case: dict[str, Any], tmp_dir: Path
) -> tuple[str, str]:
    """Run real `evaluate_candidates` + `filter_decisions_from_ranked` and build Lean `#eval` expr."""
    cid = case["case_id"]
    tsv_path = tmp_dir / f"{cid}.tsv"
    header = "smiles\tvendor\tcatalog_id\tpack_value\tpack_unit\tprice\tcurrency\tin_stock\tsource_id\tsource_url\tlicense"
    lines = [header]
    for idx, r in enumerate(case["rows"]):
        pv = "" if r["pack_val"] is None else str(r["pack_val"])
        pu = "" if r["pack_unit"] is None else str(r["pack_unit"])
        pr = "" if r["price"] is None else str(r["price"])
        cu = "" if r["currency"] is None else str(r["currency"])
        st = "" if r["in_stock"] is None else ("true" if r["in_stock"] else "false")
        lines.append(
            f"{r['smiles']}\t{r['vendor']}\tCAT-{idx}\t{pv}\t{pu}\t{pr}\t{cu}\t{st}\t"
            f"local_stock\thttps://example.org/{cid}/{idx}\tCC0-1.0"
        )
    tsv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    ranked, records = evaluate_candidates(
        [case["smiles"]],
        max_steps=int(case["max_steps"]),
        backends=("local_stock",),
        route_backend="template",
        planner="lookup_smarts",
        offline=True,
        cache_path=tmp_dir / f"{cid}_cache.sqlite",
        stock_paths=[tsv_path],
        reaction_db_path=tmp_dir / "no_rxn.sqlite",
        scscore_weights=tmp_dir / "no_sc.json.gz",
        include_scores=False,
        currency="USD",
    )
    budget: Optional[float] = case["budget"]
    summary = filter_decisions_from_ranked(
        ranked,
        budget=budget,
        budget_currency="USD",
    )
    res = ranked[0]
    rec = records[0]
    dec = summary.decisions[0]

    if dec.keep:
        py_expected = "OK"
    elif res.status not in {"buyable", "makeable"}:
        py_expected = "acyclic_atom_balanced_route"
    else:
        py_expected = "strict_homogeneous_budget"

    budget_cents = None if budget is None else int(round(budget * 100.0))
    t_vec = atom_vec_from_smiles(rec.molecule.smiles)

    if res.status == "buyable":
        is_direct_buy = True
        steps = 0
        bond_acyclic = True
        d_vec = AtomCountVec()
        r_vecs: tuple[AtomCountVec, ...] = ()
        all_stock = True
        has_price = res.best_price is not None and res.best_price.value is not None
        price_cents = int(round(float(res.best_price.value) * 100.0)) if has_price and res.best_price else None
        cost_complete = has_price
        is_usdg = bool(res.best_price is not None and res.best_price.unit == "USD/g")
    elif res.status == "makeable":
        solved_routes = [rt for rt in rec.routes if rt.leaves and all(rt.leaves_purchasable)]
        best_rt = solved_routes[0]
        is_direct_buy = False
        steps = best_rt.n_steps
        d_acc = AtomCountVec()
        bond_acyclic = True
        for st_obj in best_rt.steps:
            tpl_name = (st_obj.template_or_reference or "").split(":", 1)[0].strip()
            step_ok = ("!@" in (st_obj.template_or_reference or "")) or (
                tpl_name in UNIMOLECULAR_RING_OPENING_TEMPLATES and len(st_obj.reactants) == 1
            )
            bond_acyclic = bond_acyclic and step_ok
            d_acc = d_acc.add(atom_vec_from_counts(RETRO_TEMPLATE_DELTAS[tpl_name]))
        d_vec = d_acc
        r_vecs = tuple(atom_vec_from_smiles(lf) for lf in best_rt.leaves)
        all_stock = all(best_rt.leaves_purchasable)
        has_price = res.est_cost is not None and res.est_cost.value is not None
        price_cents = int(round(float(res.est_cost.value) * 100.0)) if has_price and res.est_cost else None
        cost_complete = bool(res.cost_complete)
        is_usdg = bool(res.est_cost is not None and res.est_cost.unit == "USD/g")
    else:
        is_direct_buy = False
        steps = 1
        bond_acyclic = True
        d_vec = AtomCountVec()
        r_vecs = (t_vec,)
        all_stock = False
        has_price = res.best_price is not None and res.best_price.value is not None
        price_cents = int(round(float(res.best_price.value) * 100.0)) if has_price and res.best_price else None
        cost_complete = False
        is_usdg = bool(res.best_price is not None and res.best_price.unit == "USD/g")

    r_list_lean = "[" + ", ".join(_atom_vec_lean_lit(v) for v in r_vecs) + "]"
    prop_lean = (
        f"{{ isDirectBuy := {_bool_lean(is_direct_buy)}, "
        f"steps := {steps}, maxSteps := {int(case['max_steps'])}, "
        f"bondIsAcyclic := {_bool_lean(bond_acyclic)}, "
        f"targetAtoms := {_atom_vec_lean_lit(t_vec)}, "
        f"leavingDelta := {_atom_vec_lean_lit(d_vec)}, "
        f"reactantAtoms := {r_list_lean}, "
        f"allLeavesInStock := {_bool_lean(all_stock)}, "
        f"costComplete := {_bool_lean(cost_complete)}, "
        f"isUsdPerGram := {_bool_lean(is_usdg)}, "
        f"priceCents := {_opt_nat_lean(price_cents)}, "
        f"budgetCents := {_opt_nat_lean(budget_cents)}, "
        f"sourceLicenses := [.permissive], "
        f"commercialMode := false }}"
    )
    lean_expr = (
        f'("{cid}", match BuyOrMakeLean.checkBuyOrMake {prop_lean} with '
        f'| none => "OK" | some v => v.rule)'
    )
    return py_expected, lean_expr


def _python_eval_moq_case(case: dict[str, Any]) -> tuple[Optional[int], bool, bool]:
    req_ug = int(case["req_ug"])
    req_b_ug = int(case["req_b_ug"])
    skus: list[tuple[int, int]] = case["skus"]
    max_outlay_cents: Optional[int] = case["max_outlay_cents"]
    require_complete: bool = case["require_complete"]

    offers: list[Offer] = []
    for idx, (pack_ug, price_cents) in enumerate(skus):
        offers.append(
            Offer(
                inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                vendor=f"V{idx}",
                pack_size=Quantity(value=float(pack_ug), unit="ug", raw_unit_str=f"{pack_ug} ug"),
                price=float(price_cents) / 100.0,
                currency="USD",
                price_basis="listed_catalog",
                source_id="local_stock",
                source_url="https://example.org/sku",
                license="CC0-1.0",
            )
        )
    outlay_q, _, _ = min_outlay_for_offers(
        offers,
        currency="USD",
        required_grams=float(req_ug) * 1e-6,
    )
    outlay_cents = int(round(outlay_q.value * 100.0)) if (outlay_q and outlay_q.value is not None) else None

    if max_outlay_cents is None:
        gate_ok = True
    elif outlay_cents is None:
        gate_ok = not require_complete
    else:
        gate_ok = outlay_cents <= max_outlay_cents

    # Check shared-synthon subadditivity on first SKU when present
    subadditive_ok = True
    if skus:
        p0, c0 = skus[0]
        pooled = math.ceil((req_ug + req_b_ug) / p0) * c0
        split = math.ceil(req_ug / p0) * c0 + math.ceil(req_b_ug / p0) * c0
        subadditive_ok = pooled <= split

    return outlay_cents, gate_ok, subadditive_ok


def _python_compound_yield_ug(base_ug: int, yield_pm: int, depth: int) -> int:
    cur = base_ug
    for _ in range(depth):
        cur = (cur * 1000 + yield_pm - 1) // yield_pm
    return cur


def run_lean_verification(*, verbose: bool = True) -> dict[str, Any]:
    lake_exe = find_lake_executable()
    if lake_exe is None:
        raise RuntimeError("Lean `lake` executable not found on PATH or ~/.elan/bin")

    # 1. Build BuyOrMakeLean
    build_proc = subprocess.run(
        [lake_exe, "build"],
        cwd=str(LEAN_DIR),
        capture_output=True,
        text=True,
        check=False,
    )
    if build_proc.returncode != 0:
        raise RuntimeError(
            f"`lake build` failed (exit {build_proc.returncode}):\n"
            f"{build_proc.stdout}\n{build_proc.stderr}"
        )

    # 2. Check `#print axioms` on every theorem in BuyOrMakeLean.lean
    theorem_names = extract_theorem_names(LEAN_FILE)
    axiom_lines = ["import BuyOrMakeLean"] + [
        f"#print axioms BuyOrMakeLean.{name}" for name in theorem_names
    ]

    # 3. Generate >=120 differential parity `#eval` expressions
    suites = build_differential_cases()
    expected_records: list[dict[str, Any]] = []
    eval_exprs: list[str] = []

    for c in suites["proposals"]:
        py_rule = _python_check_proposal(c)
        expected_records.append({"id": c["case_id"].strip(), "kind": "proposal", "expected": py_rule or "OK"})
        r_list_lean = "[" + ", ".join(_atom_vec_lean_lit(v) for v in c["reactant_atoms"]) + "]"
        lic_list_lean = "[" + ", ".join(f".{lic}" for lic in c["licenses"]) + "]"
        prop_lean = (
            f"{{ isDirectBuy := {_bool_lean(c['is_direct_buy'])}, "
            f"steps := {c['steps']}, maxSteps := {c['max_steps']}, "
            f"bondIsAcyclic := {_bool_lean(c['bond_is_acyclic'])}, "
            f"targetAtoms := {_atom_vec_lean_lit(c['target_atoms'])}, "
            f"leavingDelta := {_atom_vec_lean_lit(c['leaving_delta'])}, "
            f"reactantAtoms := {r_list_lean}, "
            f"allLeavesInStock := {_bool_lean(c['all_leaves_in_stock'])}, "
            f"costComplete := {_bool_lean(c['cost_complete'])}, "
            f"isUsdPerGram := {_bool_lean(c['is_usd_per_gram'])}, "
            f"priceCents := {_opt_nat_lean(c['price_cents'])}, "
            f"budgetCents := {_opt_nat_lean(c['budget_cents'])}, "
            f"sourceLicenses := {lic_list_lean}, "
            f"commercialMode := {_bool_lean(c['commercial_mode'])} }}"
        )
        cid = c["case_id"].strip()
        eval_exprs.append(
            f'("{cid}", match BuyOrMakeLean.checkBuyOrMake {prop_lean} with '
            f'| none => "OK" | some v => v.rule)'
        )

    for c in suites["moq"]:
        outlay_cents, gate_ok, sub_ok = _python_eval_moq_case(c)
        out_str = "NONE" if outlay_cents is None else str(outlay_cents)
        exp_str = f"{out_str}|{'T' if gate_ok else 'F'}|{'T' if sub_ok else 'F'}"
        cid = c["case_id"]
        expected_records.append({"id": cid, "kind": "moq", "expected": exp_str})
        sku_lean = "[" + ", ".join(f"({p}, {pr})" for p, pr in c["skus"]) + "]"
        p0, c0 = c["skus"][0] if c["skus"] else (1_000_000, 0)
        eval_exprs.append(
            f'("{cid}", '
            f"let outOpt := BuyOrMakeLean.minSkuOutlayCents {c['req_ug']} {sku_lean}; "
            f"let gate := BuyOrMakeLean.moqOutlayPolicyOK outOpt {_opt_nat_lean(c['max_outlay_cents'])} "
            f"{_bool_lean(c['require_complete'])}; "
            f"let sub := decide (BuyOrMakeLean.discretePackOutlayCents ({c['req_ug']} + {c['req_b_ug']}) {p0} {c0} <= "
            f"BuyOrMakeLean.discretePackOutlayCents {c['req_ug']} {p0} {c0} + "
            f"BuyOrMakeLean.discretePackOutlayCents {c['req_b_ug']} {p0} {c0}); "
            f'let outS := match outOpt with | none => "NONE" | some v => toString v; '
            f'let gS := if gate then "T" else "F"; '
            f'let sS := if sub then "T" else "F"; '
            f's!"{{outS}}|{{gS}}|{{sS}}")'
        )

    for c in suites["yield_pareto"]:
        py_mass = _python_compound_yield_ug(c["base_ug"], c["yield_pm"], c["depth"])
        ac, al, ap = c["a"]
        bc, bl, bp = c["b"]
        opt_a = GreenProcurementOption(
            option_id="a",
            is_direct_buy=True,
            steps=0,
            max_steps=1,
            bond_is_acyclic=True,
            target_atoms=AtomCountVec(),
            leaving_delta=AtomCountVec(),
            reactant_atoms=(),
            all_leaves_in_stock=True,
            cost_complete=True,
            is_usd_per_gram=True,
            usd_per_g=ac / 100.0,
            lead_time_days=al,
            pmi_kg_per_kg=ap / 1000.0,
        )
        opt_b = GreenProcurementOption(
            option_id="b",
            is_direct_buy=True,
            steps=0,
            max_steps=1,
            bond_is_acyclic=True,
            target_atoms=AtomCountVec(),
            leaving_delta=AtomCountVec(),
            reactant_atoms=(),
            all_leaves_in_stock=True,
            cost_complete=True,
            is_usd_per_gram=True,
            usd_per_g=bc / 100.0,
            lead_time_days=bl,
            pmi_kg_per_kg=bp / 1000.0,
        )
        py_dom = pareto_dominates(opt_a, opt_b)
        exp_str = f"{py_mass}|{'T' if py_dom else 'F'}"
        cid = c["case_id"]
        expected_records.append({"id": cid, "kind": "yield_pareto", "expected": exp_str})
        eval_exprs.append(
            f'("{cid}", '
            f"let m := BuyOrMakeLean.compoundYieldMassUg {c['base_ug']} {c['yield_pm']} {c['depth']}; "
            f"let d := BuyOrMakeLean.paretoDominates "
            f"{{ costCentsPerGram := {ac}, leadTimeDays := {al}, pmiPermille := {ap} }} "
            f"{{ costCentsPerGram := {bc}, leadTimeDays := {bl}, pmiPermille := {bp} }}; "
            f'let dS := if d then "T" else "F"; '
            f's!"{{m}}|{{dS}}")'
        )

    with tempfile.TemporaryDirectory(prefix="lean_shipping_gate_") as ship_tmp_raw:
        ship_tmp = Path(ship_tmp_raw)
        for c in suites["shipping_gate"]:
            py_exp, lean_expr = _eval_shipping_gate_case(c, ship_tmp)
            expected_records.append(
                {"id": c["case_id"], "kind": "shipping_gate", "expected": py_exp}
            )
            eval_exprs.append(lean_expr)

    total_cases = len(expected_records)
    assert total_cases >= 170, f"Expected >= 170 differential cases, got {total_cases}"

    harness_lines = (
        axiom_lines
        + [
            "",
            "def parityResults : List (String × String) := [",
            ",\n".join("  " + e for e in eval_exprs),
            "]",
            'def emitParity : IO Unit := do',
            '  for (cid, val) in parityResults do',
            '    IO.println s!"PARITY::{cid}::{val}"',
            "#eval emitParity",
            "",
        ]
    )

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".lean",
        prefix="lean_parity_",
        encoding="utf-8",
        delete=False,
    ) as tmp:
        tmp.write("\n".join(harness_lines))
        tmp_path = Path(tmp.name)

    try:
        proc = subprocess.run(
            [lake_exe, "env", "lean", str(tmp_path)],
            cwd=str(LEAN_DIR),
            capture_output=True,
            text=True,
            check=False,
        )
    finally:
        if tmp_path.exists():
            try:
                os.unlink(tmp_path)
            except OSError:
                pass

    if proc.returncode != 0:
        raise RuntimeError(
            f"Lean axiom/parity harness failed (exit {proc.returncode}):\n"
            f"{proc.stdout}\n{proc.stderr}"
        )

    combined_out = proc.stdout + "\n" + proc.stderr
    axioms_by_thm: dict[str, list[str]] = {}
    for line in combined_out.splitlines():
        line_s = line.strip()
        m_none = re.match(r"^'BuyOrMakeLean\.([^']+)' does not depend on any axioms$", line_s)
        if m_none:
            axioms_by_thm[m_none.group(1)] = []
            continue
        m_dep = re.match(r"^'BuyOrMakeLean\.([^']+)' depends on axioms: \[(.*)\]$", line_s)
        if m_dep:
            ax_list = [a.strip() for a in m_dep.group(2).split(",") if a.strip()]
            axioms_by_thm[m_dep.group(1)] = ax_list

    missing_thms = [t for t in theorem_names if t not in axioms_by_thm]
    if missing_thms:
        raise AssertionError(f"#print axioms output missing for theorems: {missing_thms}")

    for thm, ax_list in axioms_by_thm.items():
        disallowed = [ax for ax in ax_list if ax not in ALLOWED_AXIOMS]
        if disallowed:
            raise AssertionError(
                f"Theorem BuyOrMakeLean.{thm} depends on disallowed axioms: {disallowed} "
                f"(allowed: {sorted(ALLOWED_AXIOMS)})"
            )

    lean_results: dict[str, str] = {}
    for line in combined_out.splitlines():
        if line.startswith("PARITY::"):
            _, cid, val = line.split("::", 2)
            lean_results[cid] = val.strip()

    mismatches: list[str] = []
    for rec in expected_records:
        cid = rec["id"]
        exp = rec["expected"]
        got = lean_results.get(cid)
        if got != exp:
            mismatches.append(f"{cid} ({rec['kind']}): python={exp!r} != lean={got!r}")

    if mismatches:
        raise AssertionError(
            f"Lean/Python differential parity failed on {len(mismatches)}/{total_cases} cases:\n"
            + "\n".join(mismatches[:20])
        )

    cases_shipping = len(suites["moq"]) + len(suites["shipping_gate"])
    cases_mirror = len(suites["proposals"]) + len(suites["yield_pareto"])

    summary = {
        "theorems_checked": len(theorem_names),
        "theorem_names": theorem_names,
        "axioms_used_union": sorted({ax for axs in axioms_by_thm.values() for ax in axs}),
        "custom_axioms": 0,
        "sorry_count": 0,
        "differential_cases_total": total_cases,
        "cases_against_shipping_code": cases_shipping,
        "cases_against_mirror": cases_mirror,
        "suite_counts": {k: len(v) for k, v in suites.items()},
        "suite_oracle_map": {
            "shipping_gate": "shipping: buyormake.rank.evaluate_candidates + filter_decisions_from_ranked",
            "moq": "shipping: buyormake.cost.min_outlay_for_offers",
            "proposals": "mirror + delegated gate: community_extensions.pareto_green_triage.verify_option_guardrails",
            "yield_pareto": "mirror: community_extensions.pareto_green_triage.pareto_dominates + fixed-point yield",
        },
        "mismatches": 0,
    }
    if verbose:
        print(json.dumps(summary, indent=2))
    return summary


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check Lean 4.36.0 BuyOrMakeLean axioms and >=120-case Python differential parity."
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress JSON summary output on success.")
    args = parser.parse_args(argv)
    run_lean_verification(verbose=not args.quiet)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
