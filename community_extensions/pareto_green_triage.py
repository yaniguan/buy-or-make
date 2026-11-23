"""3D Cost-LeadTime-Greenness (PMI / E-Factor) Pareto Router (Idea #05).

Inspired by AISLE and sustainable pharmaceutical retrosynthesis routing threads.
Computes exact non-dominated 3D Pareto frontiers over
`(usd_per_g, lead_time_days, pmi_kg_per_kg)` while strictly enforcing:
1. `null`-price and incomplete-cost abstention (`unpriced_never_passes_budget`),
2. Homogeneous `USD/g` unit verification,
3. Acyclic bond cleavage and exact heavy-atom conservation (`BuyOrMakeLean`).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AtomCountVec:
    """Exact heavy-atom vector (C, N, O, S, Hal) mirroring `BuyOrMakeLean.AtomVec`."""

    c: int
    n: int
    o: int
    s: int
    hal: int

    def add(self, other: AtomCountVec) -> AtomCountVec:
        return AtomCountVec(
            c=self.c + other.c,
            n=self.n + other.n,
            o=self.o + other.o,
            s=self.s + other.s,
            hal=self.hal + other.hal,
        )


def _sum_atoms(vecs: Sequence[AtomCountVec]) -> AtomCountVec:
    acc = AtomCountVec(0, 0, 0, 0, 0)
    for vec in vecs:
        acc = acc.add(vec)
    return acc


@dataclass(frozen=True)
class GreenProcurementOption:
    """Direct-buy or N-step synthesis option annotated with cost, lead time, and PMI."""

    option_id: str
    is_direct_buy: bool
    steps: int
    max_steps: int
    bond_is_acyclic: bool
    target_atoms: AtomCountVec
    leaving_delta: AtomCountVec
    reactant_atoms: tuple[AtomCountVec, ...]
    all_leaves_in_stock: bool
    cost_complete: bool
    is_usd_per_gram: bool
    usd_per_g: float | None
    lead_time_days: int
    pmi_kg_per_kg: float

    @property
    def e_factor_kg_per_kg(self) -> float:
        return max(0.0, float(self.pmi_kg_per_kg) - 1.0)


@dataclass(frozen=True)
class ParetoTriageDecision:
    """Certified triage output separating valid Pareto options from quarantined ones."""

    selected_option: GreenProcurementOption | None
    pareto_frontier: tuple[GreenProcurementOption, ...]
    abstained_unpriced_ids: tuple[str, ...]
    rejected_stoichiometry_ids: tuple[str, ...]
    over_budget_ids: tuple[str, ...]


def verify_option_guardrails(
    opt: GreenProcurementOption,
    budget_usd_per_g: float | None = None,
) -> tuple[bool, str | None]:
    """Exact Python mirror of `BuyOrMakeLean.checkBuyOrMake` stoichiometry & budget checks."""
    if (
        opt.lead_time_days < 0
        or opt.pmi_kg_per_kg < 1.0
        or not math.isfinite(opt.pmi_kg_per_kg)
    ):
        return False, "invalid_physical_metrics"

    if opt.is_direct_buy:
        if opt.steps != 0:
            return False, "acyclic_atom_balanced_route"
    else:
        if not (1 <= opt.steps <= opt.max_steps):
            return False, "acyclic_atom_balanced_route"
        if not opt.bond_is_acyclic or not opt.all_leaves_in_stock:
            return False, "acyclic_atom_balanced_route"
        if _sum_atoms(opt.reactant_atoms) != opt.target_atoms.add(opt.leaving_delta):
            return False, "acyclic_atom_balanced_route"

    if opt.usd_per_g is None or not opt.cost_complete or not opt.is_usd_per_gram:
        return False, "strict_homogeneous_budget"
    if opt.usd_per_g < 0.0 or not math.isfinite(opt.usd_per_g):
        return False, "strict_homogeneous_budget"
    if budget_usd_per_g is not None and opt.usd_per_g > budget_usd_per_g + 1e-12:
        return False, "over_budget"

    return True, None


def pareto_dominates(a: GreenProcurementOption, b: GreenProcurementOption) -> bool:
    """Return True iff `a` is no worse than `b` on (cost, lead_time, PMI) and strictly better on >=1."""
    if a.usd_per_g is None or b.usd_per_g is None:
        raise ValueError("pareto_dominates requires priced options")
    no_worse = (
        a.usd_per_g <= b.usd_per_g + 1e-12
        and a.lead_time_days <= b.lead_time_days
        and a.pmi_kg_per_kg <= b.pmi_kg_per_kg + 1e-12
    )
    strictly_better = (
        a.usd_per_g < b.usd_per_g - 1e-12
        or a.lead_time_days < b.lead_time_days
        or a.pmi_kg_per_kg < b.pmi_kg_per_kg - 1e-12
    )
    return bool(no_worse and strictly_better)


def compute_pareto_frontier(
    options: Sequence[GreenProcurementOption],
    budget_usd_per_g: float | None = None,
) -> ParetoTriageDecision:
    """Filter stoichiometrically/budget-valid options and compute exact 3D Pareto frontier."""
    valid: list[GreenProcurementOption] = []
    abstained_unpriced: list[str] = []
    rejected_stoich: list[str] = []
    over_budget: list[str] = []

    for opt in options:
        ok, rule = verify_option_guardrails(opt, budget_usd_per_g=budget_usd_per_g)
        if ok:
            valid.append(opt)
        elif rule == "acyclic_atom_balanced_route":
            rejected_stoich.append(opt.option_id)
        elif rule == "strict_homogeneous_budget":
            abstained_unpriced.append(opt.option_id)
        elif rule == "over_budget":
            over_budget.append(opt.option_id)
        else:
            rejected_stoich.append(opt.option_id)

    frontier: list[GreenProcurementOption] = []
    for cand in valid:
        if any(
            pareto_dominates(other, cand)
            for other in valid
            if other.option_id != cand.option_id
        ):
            continue
        frontier.append(cand)

    frontier.sort(
        key=lambda o: (float(o.usd_per_g or 0.0), o.lead_time_days, o.pmi_kg_per_kg)
    )
    return ParetoTriageDecision(
        selected_option=None,
        pareto_frontier=tuple(frontier),
        abstained_unpriced_ids=tuple(abstained_unpriced),
        rejected_stoichiometry_ids=tuple(rejected_stoich),
        over_budget_ids=tuple(over_budget),
    )


def select_pareto_green_option(
    options: Sequence[GreenProcurementOption],
    budget_usd_per_g: float = 60.0,
    max_lead_days: int = 10,
) -> ParetoTriageDecision:
    """Select the balanced knee on the 3D Pareto frontier respecting budget and lead time."""
    base = compute_pareto_frontier(options, budget_usd_per_g=budget_usd_per_g)
    if not base.pareto_frontier:
        return base

    on_time = [o for o in base.pareto_frontier if o.lead_time_days <= max_lead_days]
    pool = on_time if on_time else list(base.pareto_frontier)

    # Min-max utopia-point Euclidean distance on normalized (cost, lead_time, PMI)
    costs = np.array([float(o.usd_per_g or 0.0) for o in pool], dtype=np.float64)
    leads = np.array([float(o.lead_time_days) for o in pool], dtype=np.float64)
    pmis = np.array([float(o.pmi_kg_per_kg) for o in pool], dtype=np.float64)

    c_norm = (costs - costs.min()) / max(float(costs.max() - costs.min()), 1.0)
    l_norm = (leads - leads.min()) / max(float(leads.max() - leads.min()), 1.0)
    p_norm = (pmis - pmis.min()) / max(float(pmis.max() - pmis.min()), 1.0)

    # Weight greenness (PMI) and lead-time reliability alongside unit cost
    utopia_score = 0.30 * (c_norm**2) + 0.30 * (l_norm**2) + 0.40 * (p_norm**2)
    best_idx = int(np.argmin(utopia_score))
    return ParetoTriageDecision(
        selected_option=pool[best_idx],
        pareto_frontier=base.pareto_frontier,
        abstained_unpriced_ids=base.abstained_unpriced_ids,
        rejected_stoichiometry_ids=base.rejected_stoichiometry_ids,
        over_budget_ids=base.over_budget_ids,
    )


def select_scalar_cost_only(
    options: Sequence[GreenProcurementOption],
    budget_usd_per_g: float = 60.0,
) -> GreenProcurementOption | None:
    """Baseline cost-only selector picking minimum USD/g among guardrail-valid options."""
    valid = [
        o
        for o in options
        if verify_option_guardrails(o, budget_usd_per_g=budget_usd_per_g)[0]
    ]
    if not valid:
        return None
    return min(valid, key=lambda o: (float(o.usd_per_g or 0.0), o.option_id))


def run_benchmark(
    n_candidates: int = 80,
    budget_usd_per_g: float = 60.0,
    max_lead_days: int = 10,
    seed: int = 42,
) -> dict[str, object]:
    """Compare `scalar_cost_only` vs `pareto_green_router` across 80 procurement targets."""
    if n_candidates < 1:
        raise ValueError("n_candidates must be >= 1")
    rng = np.random.default_rng(seed)

    scalar_pmis: list[float] = []
    scalar_on_time = 0
    scalar_within_budget = 0

    pareto_pmis: list[float] = []
    pareto_on_time = 0
    pareto_within_budget = 0
    null_price_leaks = 0
    stoichiometry_leaks = 0

    target_atoms = AtomCountVec(c=14, n=2, o=2, s=0, hal=0)
    leaving_hcl = AtomCountVec(c=0, n=0, o=0, s=0, hal=1)
    leaving_h2o = AtomCountVec(c=0, n=0, o=1, s=0, hal=0)
    reactants_chlorinated = (
        AtomCountVec(c=8, n=0, o=1, s=0, hal=1),
        AtomCountVec(c=6, n=2, o=1, s=0, hal=0),
    )
    reactants_green_amide = (
        AtomCountVec(c=8, n=0, o=2, s=0, hal=0),
        AtomCountVec(c=6, n=2, o=1, s=0, hal=0),
    )

    for idx in range(n_candidates):
        # Option A: Cheap stoichiometric acyl-chloride route (high chlorinated PMI, slow shipping)
        opt_cheap_dirty = GreenProcurementOption(
            option_id=f"cand_{idx}_cheap_dirty",
            is_direct_buy=False,
            steps=2,
            max_steps=3,
            bond_is_acyclic=True,
            target_atoms=target_atoms,
            leaving_delta=leaving_hcl,
            reactant_atoms=reactants_chlorinated,
            all_leaves_in_stock=True,
            cost_complete=True,
            is_usd_per_gram=True,
            usd_per_g=round(float(rng.uniform(14.0, 24.0)), 2),
            lead_time_days=int(rng.integers(11, 19)),
            pmi_kg_per_kg=round(float(rng.uniform(68.0, 112.0)), 2),
        )
        # Option B: Catalytic aqueous/bio-solvent 1-step route (moderate cost, fast, low PMI)
        opt_green_fast = GreenProcurementOption(
            option_id=f"cand_{idx}_green_catalytic",
            is_direct_buy=False,
            steps=1,
            max_steps=3,
            bond_is_acyclic=True,
            target_atoms=target_atoms,
            leaving_delta=leaving_h2o,
            reactant_atoms=reactants_green_amide,
            all_leaves_in_stock=True,
            cost_complete=True,
            is_usd_per_gram=True,
            usd_per_g=round(float(rng.uniform(26.0, 44.0)), 2),
            lead_time_days=int(rng.integers(3, 8)),
            pmi_kg_per_kg=round(float(rng.uniform(13.0, 26.0)), 2),
        )
        # Option C: Direct domestic catalog stock (higher price, 2-day delivery, low repackaging PMI)
        opt_direct_stock = GreenProcurementOption(
            option_id=f"cand_{idx}_direct_buy",
            is_direct_buy=True,
            steps=0,
            max_steps=3,
            bond_is_acyclic=True,
            target_atoms=target_atoms,
            leaving_delta=AtomCountVec(0, 0, 0, 0, 0),
            reactant_atoms=(),
            all_leaves_in_stock=True,
            cost_complete=True,
            is_usd_per_gram=True,
            usd_per_g=round(float(rng.uniform(46.0, 58.0)), 2),
            lead_time_days=2,
            pmi_kg_per_kg=round(float(rng.uniform(18.0, 30.0)), 2),
        )
        # Option D: Unpriced public vendor listing (`usd_per_g=None` -> must abstain!)
        opt_unpriced = GreenProcurementOption(
            option_id=f"cand_{idx}_unpriced",
            is_direct_buy=True,
            steps=0,
            max_steps=3,
            bond_is_acyclic=True,
            target_atoms=target_atoms,
            leaving_delta=AtomCountVec(0, 0, 0, 0, 0),
            reactant_atoms=(),
            all_leaves_in_stock=True,
            cost_complete=False,
            is_usd_per_gram=False,
            usd_per_g=None,
            lead_time_days=1,
            pmi_kg_per_kg=5.0,
        )
        # Option E: Endocyclic lactam ring-cleavage trap (`bond_is_acyclic=False` -> must reject!)
        opt_ring_trap = GreenProcurementOption(
            option_id=f"cand_{idx}_ring_trap",
            is_direct_buy=False,
            steps=1,
            max_steps=3,
            bond_is_acyclic=False,
            target_atoms=target_atoms,
            leaving_delta=leaving_h2o,
            reactant_atoms=reactants_green_amide,
            all_leaves_in_stock=True,
            cost_complete=True,
            is_usd_per_gram=True,
            usd_per_g=8.0,
            lead_time_days=2,
            pmi_kg_per_kg=8.0,
        )

        pool = (
            opt_cheap_dirty,
            opt_green_fast,
            opt_direct_stock,
            opt_unpriced,
            opt_ring_trap,
        )

        s_pick = select_scalar_cost_only(pool, budget_usd_per_g=budget_usd_per_g)
        assert s_pick is not None and s_pick.usd_per_g is not None
        scalar_pmis.append(s_pick.pmi_kg_per_kg)
        if s_pick.lead_time_days <= max_lead_days:
            scalar_on_time += 1
        if s_pick.usd_per_g <= budget_usd_per_g:
            scalar_within_budget += 1

        p_res = select_pareto_green_option(
            pool, budget_usd_per_g=budget_usd_per_g, max_lead_days=max_lead_days
        )
        p_pick = p_res.selected_option
        assert p_pick is not None
        if p_pick.usd_per_g is None or not p_pick.cost_complete:
            null_price_leaks += 1
        if not p_pick.is_direct_buy and not p_pick.bond_is_acyclic:
            stoichiometry_leaks += 1

        pareto_pmis.append(p_pick.pmi_kg_per_kg)
        if p_pick.lead_time_days <= max_lead_days:
            pareto_on_time += 1
        if p_pick.usd_per_g is not None and p_pick.usd_per_g <= budget_usd_per_g:
            pareto_within_budget += 1

    mean_scalar_pmi = float(np.mean(scalar_pmis))
    mean_pareto_pmi = float(np.mean(pareto_pmis))
    waste_reduction_pct = (
        100.0 * (mean_scalar_pmi - mean_pareto_pmi) / max(mean_scalar_pmi, 1e-9)
    )

    return {
        "idea_id": "idea_05_pareto_green_triage",
        "n_candidates": n_candidates,
        "budget_usd_per_g": budget_usd_per_g,
        "max_lead_days": max_lead_days,
        "scalar_cost_only": {
            "mean_pmi_kg_per_kg": round(mean_scalar_pmi, 2),
            "on_time_delivery_rate": round(scalar_on_time / n_candidates, 4),
            "within_budget_rate": round(scalar_within_budget / n_candidates, 4),
        },
        "pareto_green_router": {
            "mean_pmi_kg_per_kg": round(mean_pareto_pmi, 2),
            "mean_pmi_waste_reduction_pct": round(waste_reduction_pct, 2),
            "on_time_delivery_rate": round(pareto_on_time / n_candidates, 4),
            "within_budget_rate": round(pareto_within_budget / n_candidates, 4),
            "null_price_leaks": null_price_leaks,
            "stoichiometry_leaks": stoichiometry_leaks,
        },
    }
