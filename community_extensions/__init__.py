"""Community research extensions for buy-or-make (L5 procurement & retrosynthesis triage).

Modules:
- pareto_green_triage: 3D Cost-LeadTime-Greenness (PMI / E-Factor) Pareto Router (Idea #05).
- plate_pack_procurement: Shared-Synthon Combinatorial Cart Optimizer (Idea #06).
"""

from __future__ import annotations

from .pareto_green_triage import (
    AtomCountVec,
    GreenProcurementOption,
    ParetoTriageDecision,
    atom_vec_from_counts,
    atom_vec_from_mol,
    atom_vec_from_smiles,
    compute_pareto_frontier,
    pareto_dominates,
    select_pareto_green_option,
    select_scalar_cost_only,
    verify_option_guardrails,
)
from .pareto_green_triage import (
    run_benchmark as run_pareto_green_benchmark,
)
from .plate_pack_procurement import (
    CartOptimizationResult,
    SynthonPack,
    TargetAnalogueSpec,
    optimize_independent_per_target_cart,
    optimize_shared_synthon_cart,
)
from .plate_pack_procurement import (
    run_benchmark as run_plate_pack_benchmark,
)

__all__ = [
    "AtomCountVec",
    "CartOptimizationResult",
    "GreenProcurementOption",
    "ParetoTriageDecision",
    "SynthonPack",
    "TargetAnalogueSpec",
    "atom_vec_from_counts",
    "atom_vec_from_mol",
    "atom_vec_from_smiles",
    "compute_pareto_frontier",
    "optimize_independent_per_target_cart",
    "optimize_shared_synthon_cart",
    "pareto_dominates",
    "run_pareto_green_benchmark",
    "run_plate_pack_benchmark",
    "select_pareto_green_option",
    "select_scalar_cost_only",
    "verify_option_guardrails",
]
