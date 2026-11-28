"""3D Cost-LeadTime-Greenness (PMI / E-Factor) Pareto Router (Idea #05).

Inspired by AISLE and sustainable pharmaceutical retrosynthesis routing threads.
Computes exact non-dominated 3D Pareto frontiers over
`(usd_per_g, lead_time_days, pmi_kg_per_kg)` while strictly enforcing:
1. `null`-price and incomplete-cost abstention (`unpriced_never_passes_budget`),
2. Homogeneous `USD/g` unit verification,
3. Acyclic bond cleavage and 13-channel heavy-atom conservation (`BuyOrMakeLean.AtomVec`).
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

from buyormake.rank import filter_decisions_from_ranked
from buyormake.schema import Quantity, RankResult

_METAL_SYMBOLS: frozenset[str] = frozenset(
    {
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",
        "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
        "Au", "Hg", "Tl", "Pb", "Bi", "Fr", "Ra", "Ac", "Th", "Pa", "U",
    }
)


@dataclass(frozen=True)
class AtomCountVec:
    """Exact 13-channel heavy-atom vector mirroring `BuyOrMakeLean.AtomVec`.

    Separates `F`, `Cl`, `Br`, `I` individually and tracks `B`, `P`, `Si`, `metal`,
    and `other_heavy` (plus exact `extra_elements` identity when constructed from SMILES)
    so cross-halogen swaps and missing Suzuki boron / organometallic synthons are caught.
    """

    c: int = 0
    n: int = 0
    o: int = 0
    s: int = 0
    f: int = 0
    cl: int = 0
    br: int = 0
    i: int = 0
    b: int = 0
    p: int = 0
    si: int = 0
    metal: int = 0
    other_heavy: int = 0
    extra_elements: tuple[tuple[str, int], ...] = ()

    @property
    def hal(self) -> int:
        """Total halogen count (`F + Cl + Br + I`) for convenience inspection."""
        return self.f + self.cl + self.br + self.i

    def add(self, other: AtomCountVec) -> AtomCountVec:
        merged_extra: Counter[str] = Counter(dict(self.extra_elements))
        for sym, cnt in other.extra_elements:
            merged_extra[sym] += cnt
        extra_tuple = tuple(sorted((k, int(v)) for k, v in merged_extra.items() if v != 0))
        return AtomCountVec(
            c=self.c + other.c,
            n=self.n + other.n,
            o=self.o + other.o,
            s=self.s + other.s,
            f=self.f + other.f,
            cl=self.cl + other.cl,
            br=self.br + other.br,
            i=self.i + other.i,
            b=self.b + other.b,
            p=self.p + other.p,
            si=self.si + other.si,
            metal=self.metal + other.metal,
            other_heavy=self.other_heavy + other.other_heavy,
            extra_elements=extra_tuple,
        )


def atom_vec_from_counts(counts: Mapping[str, int]) -> AtomCountVec:
    """Project an element-symbol mapping (e.g., `RETRO_TEMPLATE_DELTAS['suzuki_biaryl']`) to `AtomCountVec`."""
    c = n = o = s = f = cl = br = i = b = p = si = metal = other_heavy = 0
    extra: Counter[str] = Counter()
    for raw_sym, raw_cnt in counts.items():
        cnt = int(raw_cnt)
        if cnt == 0:
            continue
        sym = raw_sym.strip()
        if sym == "H":
            continue
        if sym == "C":
            c += cnt
        elif sym == "N":
            n += cnt
        elif sym == "O":
            o += cnt
        elif sym == "S":
            s += cnt
        elif sym == "F":
            f += cnt
        elif sym == "Cl":
            cl += cnt
        elif sym == "Br":
            br += cnt
        elif sym == "I":
            i += cnt
        elif sym == "B":
            b += cnt
        elif sym == "P":
            p += cnt
        elif sym == "Si":
            si += cnt
        elif sym in _METAL_SYMBOLS:
            metal += cnt
            extra[sym] += cnt
        else:
            other_heavy += cnt
            extra[sym] += cnt
    extra_tuple = tuple(sorted((k, int(v)) for k, v in extra.items() if v != 0))
    return AtomCountVec(
        c=c,
        n=n,
        o=o,
        s=s,
        f=f,
        cl=cl,
        br=br,
        i=i,
        b=b,
        p=p,
        si=si,
        metal=metal,
        other_heavy=other_heavy,
        extra_elements=extra_tuple,
    )


def atom_vec_from_mol(mol: Chem.Mol) -> AtomCountVec:
    """Convert an RDKit `Chem.Mol` into an exact heavy-atom `AtomCountVec`."""
    counts = Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    return atom_vec_from_counts(counts)


def atom_vec_from_smiles(smiles: str) -> AtomCountVec:
    """Parse SMILES via RDKit and return its exact 13-channel `AtomCountVec`."""
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("smiles must be a non-empty string")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    return atom_vec_from_mol(mol)


def _sum_atoms(vecs: Sequence[AtomCountVec]) -> AtomCountVec:
    acc = AtomCountVec()
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
    """Verify `BuyOrMakeLean.checkBuyOrMake` stoichiometry and delegate budget gating to `buyormake.rank`."""
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

    if opt.usd_per_g is not None and (
        not math.isfinite(opt.usd_per_g) or opt.usd_per_g < 0.0
    ):
        return False, "strict_homogeneous_budget"

    qty = (
        Quantity(
            value=float(opt.usd_per_g),
            unit="USD/g" if opt.is_usd_per_gram else "USD",
        )
        if opt.usd_per_g is not None
        else None
    )
    if opt.is_direct_buy:
        rr = RankResult(
            candidate=opt.option_id,
            status="buyable",
            n_steps=0,
            best_price=qty if opt.cost_complete else None,
            cost_complete=False,
        )
    else:
        rr = RankResult(
            candidate=opt.option_id,
            status="makeable",
            n_steps=opt.steps,
            est_cost=qty,
            cost_complete=opt.cost_complete,
        )

    eff_budget = (
        float(budget_usd_per_g)
        if budget_usd_per_g is not None
        else max(0.0, float(opt.usd_per_g or 0.0))
    )
    dec = filter_decisions_from_ranked(
        [rr],
        budget=eff_budget,
        budget_currency="USD",
        require_complete_cost=True,
    ).decisions[0]
    if dec.keep:
        return True, None
    if dec.over_budget:
        return False, "over_budget"
    return False, "strict_homogeneous_budget"


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

    costs = np.array([float(o.usd_per_g or 0.0) for o in pool], dtype=np.float64)
    leads = np.array([float(o.lead_time_days) for o in pool], dtype=np.float64)
    pmis = np.array([float(o.pmi_kg_per_kg) for o in pool], dtype=np.float64)

    c_norm = (costs - costs.min()) / max(float(costs.max() - costs.min()), 1.0)
    l_norm = (leads - leads.min()) / max(float(leads.max() - leads.min()), 1.0)
    p_norm = (pmis - pmis.min()) / max(float(pmis.max() - pmis.min()), 1.0)

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

    target_atoms = AtomCountVec(c=14, n=2, o=2)
    leaving_hcl = AtomCountVec(cl=1)
    leaving_h2o = AtomCountVec(o=1)
    reactants_chlorinated = (
        AtomCountVec(c=8, o=1, cl=1),
        AtomCountVec(c=6, n=2, o=1),
    )
    reactants_green_amide = (
        AtomCountVec(c=8, o=2),
        AtomCountVec(c=6, n=2, o=1),
    )

    for idx in range(n_candidates):
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
        opt_direct_stock = GreenProcurementOption(
            option_id=f"cand_{idx}_direct_buy",
            is_direct_buy=True,
            steps=0,
            max_steps=3,
            bond_is_acyclic=True,
            target_atoms=target_atoms,
            leaving_delta=AtomCountVec(),
            reactant_atoms=(),
            all_leaves_in_stock=True,
            cost_complete=True,
            is_usd_per_gram=True,
            usd_per_g=round(float(rng.uniform(46.0, 58.0)), 2),
            lead_time_days=2,
            pmi_kg_per_kg=round(float(rng.uniform(18.0, 30.0)), 2),
        )
        opt_unpriced = GreenProcurementOption(
            option_id=f"cand_{idx}_unpriced",
            is_direct_buy=True,
            steps=0,
            max_steps=3,
            bond_is_acyclic=True,
            target_atoms=target_atoms,
            leaving_delta=AtomCountVec(),
            reactant_atoms=(),
            all_leaves_in_stock=True,
            cost_complete=False,
            is_usd_per_gram=False,
            usd_per_g=None,
            lead_time_days=1,
            pmi_kg_per_kg=5.0,
        )
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
