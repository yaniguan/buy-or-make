# `BuyOrMakeLean` -- Formal L5 Stoichiometric, MOQ/Yield & Green Pareto Guardrails in Lean 4

Self-contained Lean 4 (`leanprover/lean4:v4.36.0`, zero external dependencies, `0 sorry`, no custom axioms — `#print axioms` audited in CI to allow only standard kernel axioms `propext`, `Classical.choice`, and `Quot.sound`) formalization of the `L5` molecular procurement, retrosynthesis stoichiometry, discrete pack-size MOQ outlay, and step-depth yield compounding invariants for `buy-or-make`, synchronized with the central 15-repository formal specification in [`yaniguan/auto-lab/tree/main/lean`](https://github.com/yaniguan/auto-lab/tree/main/lean).

## Proved Theorems & Guardrails

| Component | Definition / Theorem | Guarantee |
| :--- | :--- | :--- |
| Procurement Gate | `checkBuyOrMake_ok_iff` | Exact soundness and completeness of `checkBuyOrMake` wrt `MakeOrBuyOK` |
| Acyclic Bond Rule | `ring_cleavage_rejected` | Endocyclic ring disconnections (`bondIsAcyclic = false`) never pass |
| Atom Conservation | `atom_imbalance_rejected` | Any mismatch in 13-channel `sumAtoms reactantAtoms = targetAtoms + leavingDelta` is rejected |
| Halogen Identity | `cross_halogen_swap_rejected` | Individual `F, Cl, Br, I` swaps are rejected even when lumped halogen counts match |
| Boron Stoichiometry | `suzuki_boron_imbalance_rejected` | Missing or unbalanced `B` in Suzuki-Miyaura couplings is rejected |
| Main-Group / Metals | `main_group_or_metal_imbalance_rejected` | Unbalanced `P, Si, metal, otherHeavy` synthons are rejected |
| Strict-Null Budget | `unpriced_never_passes_budget` | Missing (`none`) or incomplete leaf costs never pass an active budget gate |
| Unit Homogeneity | `non_usd_per_gram_rejected_under_budget` | Non-`USD/g` pack quotes are never mixed into per-gram budget checks |
| 3D Green Pareto | `paretoDominates_irrefl` | Irreflexive non-dominated ordering over `(cost, leadTime, PMI)` |
| MOQ Pack Coverage | `discrete_pack_covers_required_mass` | Ceiling pack count `discretePackCount reqUg packUg` always covers `reqUg` |
| Shared-Synthon MOQ | `shared_synthon_pack_count_subadditive`, `shared_synthon_outlay_subadditive` | Pooling shared building blocks across multi-branch `Retro*` trees never exceeds per-branch outlays |
| Step-Depth Yield | `yield_compounded_ge_lower_bound`, `yield_compounded_monotone_depth` | Yield-compounded mass is `>=` 100%-yield lower bound and monotone in reaction step depth |
| MOQ Outlay Gate | `moq_outlay_gate_sound` | Missing or over-ceiling discrete MOQ outlays fail strict cash-outlay gating |

## Quickstart & Differential Parity Verification

```bash
cd lean
lake build
cd ..
python scripts/check_lean_axioms_and_parity.py
```
