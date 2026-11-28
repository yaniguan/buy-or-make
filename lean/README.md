# `BuyOrMakeLean` -- Formal L5 Stoichiometric, MOQ/Yield & Green Pareto Guardrails in Lean 4

Self-contained Lean 4 (`leanprover/lean4:v4.36.0`, zero external dependencies, `0 sorry`, no custom axioms — `#print axioms` audited in CI to allow only standard kernel axioms `propext`, `Classical.choice`, and `Quot.sound`) formalization of the `L5` molecular procurement, retrosynthesis stoichiometry, discrete pack-size MOQ outlay, and step-depth yield compounding invariants for `buy-or-make`, synchronized with the central 15-repository formal specification in [`yaniguan/auto-lab/tree/main/lean`](https://github.com/yaniguan/auto-lab/tree/main/lean).

## Proved Theorems, Guardrails & Python Entry-Point Mappings

`scripts/check_lean_axioms_and_parity.py` executes **175 `#eval` differential parity cases** (`85` directly against shipping `buyormake.*` package functions + `90` against `community_extensions.pareto_green_triage`, whose budget gate also delegates to `buyormake.rank.filter_decisions_from_ranked`):

| Component | Lean Definition / Theorem | Guarantee | Python Entry Point Checked Against |
| :--- | :--- | :--- | :--- |
| Procurement Gate | `checkBuyOrMake`, `checkBuyOrMake_ok_iff` | Exact soundness and completeness of `checkBuyOrMake` wrt `MakeOrBuyOK` | `buyormake.rank.evaluate_candidates` + `buyormake.rank.filter_decisions_from_ranked` (`shipping_gate`, 40 cases) and `community_extensions.pareto_green_triage.verify_option_guardrails` (`proposals`, 60 cases) |
| Acyclic Bond Rule | `routeStoichiometryOK`, `ring_cleavage_rejected` | Endocyclic ring disconnections (`bondIsAcyclic = false`) and out-of-stock leaves (`allLeavesInStock = false`) never pass | `buyormake.routes.TemplateFallbackBackend.one_step`, `buyormake.cost.purchasable_offers`, `autolab_wire._route_is_acyclic_and_atom_balanced` |
| Atom Conservation | `atom_imbalance_rejected` | Any mismatch in 13-channel `sumAtoms reactantAtoms = targetAtoms + leavingDelta` is rejected | `buyormake.routes.TemplateFallbackBackend.one_step`, `buyormake.retro_star.RetroStarPlanner.expand_extended_smarts`, `verify_option_guardrails` |
| Halogen Identity | `cross_halogen_swap_rejected` | Individual `F, Cl, Br, I` swaps are rejected even when lumped halogen counts match | `community_extensions.pareto_green_triage.atom_vec_from_smiles` + `TemplateFallbackBackend.one_step` |
| Boron Stoichiometry | `suzuki_boron_imbalance_rejected` | Missing or unbalanced `B` in Suzuki-Miyaura couplings is rejected | `buyormake.routes.RETRO_TEMPLATE_DELTAS["suzuki_biaryl"]` + `TemplateFallbackBackend.one_step` |
| Main-Group / Metals | `main_group_or_metal_imbalance_rejected` | Unbalanced `P, Si, metal, otherHeavy` synthons are rejected | `community_extensions.pareto_green_triage.atom_vec_from_counts` / `atom_vec_from_smiles` |
| Strict-Null Budget | `budgetPolicyOK`, `unpriced_never_passes_budget` | Missing (`none`) or incomplete leaf costs never pass an active budget gate | `buyormake.rank.filter_decisions_from_ranked` & `buyormake.rank.filter_proposals` (`require_complete_cost` defaults to `True` whenever `budget` is set; explicit `require_complete_cost=False` is an opt-out) |
| Unit Homogeneity | `non_usd_per_gram_rejected_under_budget` | Non-`USD/g` quotes (raw pack `USD`, `EUR/g`, `JPY/g`) are never mixed into `USD/g` budget checks | `buyormake.cost.estimate_route_cost_bounds` (`unit="USD/g"`) + `buyormake.rank.filter_decisions_from_ranked` (`buyable_allowed_units` / `makeable_allowed_units`) |
| Transitive License | `licensePolicyOK` | `evalOnly` blocked in all modes; non-permissive blocked in `commercialMode` | `buyormake.rank._collect_provenance_and_licenses`, `buyormake.schema.select_record_license`, `proposals` parity suite |
| 3D Green Pareto | `paretoDominates`, `paretoDominates_irrefl` | Irreflexive non-dominated ordering over `(cost, leadTime, PMI)` | `community_extensions.pareto_green_triage.pareto_dominates` (`yield_pareto`, 30 cases) |
| MOQ Pack Coverage | `discretePackCount`, `discrete_pack_covers_required_mass` | Ceiling pack count `discretePackCount reqUg packUg` always covers `reqUg` | `buyormake.cost.min_outlay_for_offers` (`moq`, 45 cases) |
| Shared-Synthon MOQ | `shared_synthon_pack_count_subadditive`, `shared_synthon_outlay_subadditive` | Pooling shared building blocks across multi-branch `Retro*` trees never exceeds per-branch outlays | `buyormake.cost.min_outlay_for_offers` & `buyormake.retro_star.optimize_route_moq_outlay` (`moq`, 45 cases) |
| Step-Depth Yield | `compoundYieldMassUg`, `yield_compounded_ge_lower_bound`, `yield_compounded_monotone_depth` | Yield-compounded mass is `>=` 100%-yield lower bound and monotone in reaction step depth | `buyormake.cost.compute_leaf_depths` & `buyormake.cost.estimate_route_cost_bounds` (`yield_pareto`, 30 cases) |
| MOQ Outlay Gate | `moqOutlayPolicyOK`, `moq_outlay_gate_sound` | Missing or over-ceiling discrete MOQ outlays fail strict cash-outlay gating | `buyormake.cost.min_outlay_for_offers` + `buyormake.rank.filter_decisions_from_ranked` (`moq`, 45 cases) |

## Quickstart & Differential Parity Verification

```bash
cd lean
lake build
cd ..
python scripts/check_lean_axioms_and_parity.py
```
