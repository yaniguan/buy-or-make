# Architecture and History — `buy-or-make`

> **Auto-Lab Ecosystem Role**: Layer `L5 · Executive Planning, Procurement & Scheduling` in [Robot AI Auto-Lab (`auto-lab`)](https://github.com/yaniguan/auto-lab) ([Project Studio](https://yaniguan.github.io/auto-lab/)).

---

## 1. Object

### Problem Statement
Generative chemistry models, active-learning loops, and LLM co-scientists frequently propose candidate molecules without checking whether they can actually be bought off the shelf or synthesized in a small number of steps from commercially available building blocks—or worse, rely on hallucinated prices and opaque synthetic accessibility heuristics as proxies for real procurement cost.

`buy-or-make` solves **cost-aware molecular procurement and retrosynthetic triage**: given a list of candidate SMILES strings (including pharmaceutical salts/hydrates), it normalizes chemical structures (`RDKit`), queries local stock catalogs and live vendor registries (`PubChem PUG-REST / PUG-View`), searches for short retrosynthetic disconnections (`1`–`N` steps) back to purchasable leaves when a molecule is not directly in stock, computes conservative leaf-only stoichiometric cost lower bounds strictly from real `listed_catalog` prices, and ranks or filters proposals under explicit unit/currency-checked budget constraints while keeping missing prices `null` and unverified compounds `"unknown"`.

### Role in the Auto-Lab `L0`–`L6` Hierarchy
Within the [Auto-Lab](https://github.com/yaniguan/auto-lab) Brain/Body architecture (`L0`–`L6`), `buy-or-make` operates in **`L5 · Executive Planning, Procurement & Scheduling`**:
- **Upstream (`L6` / `L5` / `L3` Brain)**: Filters candidate molecule batches proposed by hypothesis generators ([`lab-ask-why`](https://github.com/yaniguan/lab-ask-why)), trust-region Bayesian optimizers ([`lab-decide`](https://github.com/yaniguan/lab-decide)), or molecular extraction/design workflows ([`ouroboros-ocsr`](https://github.com/yaniguan/ouroboros-ocsr), [`YLM`](https://github.com/yaniguan/YLM), [`WignerFlow`](https://github.com/yaniguan/WignerFlow)).
- **Peer (`L5`)**: Works alongside [`labshift`](https://github.com/yaniguan/labshift) via `buyormake.handoff.to_procurement_constraint` / `build_labshift_release_map` (mapping vendor `lead_time_days` and route `n_steps` to `labshift` job `release` offsets) and [`lab-automation-atlas`](https://github.com/yaniguan/lab-automation-atlas).
- **Downstream (`L0` Ledger / `L4` Verification)**: Exports Pydantic v2 `AvailabilityRecord` objects and sealed `AssayfoldRecord` envelopes (`AvailabilityRecord.to_assayfold_record()`) aligned with [`assayfold`](https://github.com/yaniguan/assayfold) (`L0`), propagating all target, route-step (`CC-BY-SA-4.0` ORD, `CC-BY-4.0` AiZynth), and leaf-offer licenses.

### Inputs, Outputs, and Scope Boundaries
- **Inputs**:
  - Candidate SMILES strings and optional molecule names (via Python API `cost_aware_rank` / `evaluate_candidates` / `filter_proposals` or CLI `buy-or-make rank|filter <file.smi>`).
  - Optional local stock files (`.tsv`, `.csv`, `.jsonl`, `.smi`, `.smi.gz`, `.sqlite`, `.hdf5`), reaction lookup SQLite index (`data/reactions.sqlite`), SCScore weights (`data/scscore/`), AiZynthFinder config (`data/aizynth/config.yml`), step cap `max_steps`, and optional per-gram `budget` (`budget_currency`, `require_complete_cost`).
- **Outputs**:
  - Sorted `list[RankResult]` classifying each candidate into `buyable | makeable | unknown | invalid`, with `name`, `n_steps`, `best_price` (`Quantity | null`), `est_cost` (`Quantity | null`), `cost_complete`, `sa_scores` (`SAScores`), human-readable `reasons`, and vendor `evidence` URLs.
  - Input-ordered `list[AvailabilityRecord]` sealed with canonical SHA-256 `record_hash`, convertible via `.to_assayfold_record()` and `to_procurement_constraint()`.
  - `FilterSummary` (`kept`, `over_budget`, `rejected`, per-candidate `ProposalDecision` with `price_status`) for co-scientist budget gating.
  - Five strict JSON Schema files (`additionalProperties: false`, zero-drift verified in CI) under `schema/`.
- **Scope Boundaries**:
  - **Not purchasing quotes**: Catalog prices vary by pack size, purity, region, currency, and institutional contracts; free public sources (PubChem, ZINC, eMolecules free, Mcule free, Enamine Open Structure Data) publish zero prices, so `best_price` and `est_cost` remain `null` unless local `listed_catalog` prices are supplied.
  - **Lower-bound stoichiometry only**: Route `est_cost` sums $1\text{ equiv}$ leaf building-block mass per $1\text{ g}$ target and excludes solvents, catalysts, reagents, non-unity yields, purification, QC, and labor.

---

## 2. Data

### Evaluated Source Registry & License Gate (`buyormake/sources.py`, `buyormake/open_corpora.py`)
Nineteen external chemical vendor, building-block, natural-product, patent-compound, reaction, and retrosynthesis benchmark sources were audited (`12` initial sources on `2026-11-18` and `7` 2024–2026 open chemistry corpora on `2026-11-25`) and encoded as validated `Source` objects (`buy-or-make sources`, with explicit `spdx_license`, ISO `YYYY-MM-DD` `license_checked`, and `license_check_method="direct"|"indirect"`):

| Source ID | Modality / Role | SPDX / License | Commercial OK | Redistribution OK | Status |
|---|---|---|---|---|---|
| `pubchem_pug` | Live vendor names, RegistryIDs & URLs (`<= 5 req/s`; no prices) | `LicenseRef-Custom-Terms` (`custom/terms`) | Yes | **No** (local SQLite cache only) | `integrated` |
| `rdkit_sascore` | Fragment + complexity heuristic (`1`–`10`) | `BSD-3-Clause` | Yes | Yes | `integrated` |
| `scscore` | Pretrained neural synthetic complexity (`1`–`5`, `6.18 MB`) | `MIT` | Yes | Yes | `integrated` |
| `zinc` | Purchasable SMILES tranches (`scripts/fetch_zinc.py`, `indirect` check) | `LicenseRef-Custom-Terms` (`custom/terms`) | Yes | **No** | `fetch_only` |
| `enamine_bb` | Building block SDFs after registration (`scripts/fetch_enamine.py`) | `LicenseRef-Custom-Terms` (`custom/terms`) | Yes (own use) | **No** | `fetch_only` |
| `emolecules` | Free `version.smi.gz` (`~350 MB`, no prices; `scripts/fetch_emolecules.py`) | `LicenseRef-Custom-Terms` (`custom/terms`) | **No** (requires ePlus) | **No** | `fetch_only` |
| `mcule` | Free in-stock building blocks (`~151 MB`; `scripts/fetch_mcule.py`) | `LicenseRef-Custom-Terms` (`custom/terms`) | Yes (own use) | **No** | `fetch_only` |
| `aizynthfinder_models` | MCTS retrosynthesis + USPTO ONNX models (`scripts/fetch_aizynth.py`) | `CC-BY-4.0` (models) / `MIT` (code) | Yes | Yes | `fetch_only` |
| `aizynthfinder_zinc_stock` | `zinc_stock_17_04_20.hdf5` (`663 MB`, ZINC-derived) | `LicenseRef-Custom-Terms` (`custom/terms`) | **No** | **No** | `fetch_only` |
| `uspto_lowe` | Lowe 1976–Sep 2016 patent reactions (`75.2 MB`; `scripts/fetch_uspto.py`) | `CC0-1.0` | Yes | Yes | `fetch_only` |
| `ord_data` | Open Reaction Database Parquet shards (`scripts/fetch_ord.py`) | `CC-BY-SA-4.0` | Yes | Yes (ShareAlike) | `fetch_only` |
| `rascore` | Reymond Group retrosynthetic classifier (`52.4 MB`) | `MIT` | Yes | Yes | `not_used` |
| `paroutes` | Multi-step patent routes (`n1`/`n5`), stopping stocks & reactions (Zenodo `7341155`) | `CC-BY-4.0` | Yes | Yes | `fetch_only` |
| `chembl` | `chembl_37` bioactive canonical SMILES & InChIKeys (May 2026, SHA-256 verified) | `CC-BY-SA-3.0` | Yes | Yes (ShareAlike) | `fetch_only` |
| `surechembl` | Biweekly patent compounds Parquet (`bulk_data/latest`, `2026-11-03`) | `CC-BY-4.0` | Yes | Yes | `fetch_only` |
| `coconut` | COCONUT 2.0 ~700k open natural products (`2026-11` S3 / Zenodo `13897048`) | `CC0-1.0` | Yes | Yes | `fetch_only` |
| `orderly` | Deduplicated, standardized ORD forward/retro Parquet benchmarks (JCIM 2024) | `CC-BY-4.0` | Yes | Yes | `fetch_only` |
| `syntheseus_uspto_50k` | Microsoft `Syntheseus` `USPTO-50K` & `Retro*` stock inventory (2024–2026) | `MIT` | Yes | Yes | `fetch_only` |
| `synplanner_data` | Standardized USPTO reactions, building blocks & SAscore targets (HF Hub) | `MIT` | Yes | Yes | `fetch_only` |

All fetched files stream into `<dest>.part`, verify pinned SHA-256 / MD5 / LFS digests via `buyormake.canonical.atomic_download`, and atomically `os.replace` into gitignored `data/` paths (`.gitignore`). Non-vendor reference corpora (`chembl`, `surechembl`, `coconut`) land in `data/reference_compounds.sqlite` (and sealed `AssayfoldRecord` envelopes) rather than `data/stock.sqlite`, preventing unpurchasable reference compounds from being mistaken for commercial vendor stock.

### Unified Pydantic v2 & JSON Schemas (`buyormake/schema.py`, `schema/`)
All domain models inherit `_Model(BaseModel)` with `ConfigDict(extra="forbid")`, enforcing finite floats (`allow_nan=False`), non-negative prices/lead times, and `price`/`currency`/`price_basis` consistency. Canonical hashing (`AvailabilityRecord.content_hash()` and `AssayfoldRecord.content_hash()`) strips `record_hash` and volatile `provenance.recorded_utc` before computing SHA-256.

### 20-Molecule Benchmark Demo (`examples/demo.py`, `examples/demo_output.md`)
Evaluates 20 curated inputs (`6` common drugs, `6` building blocks/reagents, `6` novel analogs, `2` malformed SMILES) on CPU (`AuthenticAMD`, Python `3.11.15`, RDKit `2026.03.4`):
- **Cold live PubChem run**: `38.737 s` (`1.937 s / molecule`, throttled on every `GET` and `POST`).
- **Warm offline SQLite cache run**: `0.058 s` (`0.003 s / molecule`, `~668x` speedup).
- **Triage breakdown**: `14` buyable (`70.0%`), `3` makeable (`15.0%`, including `Analog-6` via atom-conserved exocyclic Buchwald/SnAr disconnection into `2-bromopyrimidine` + `2,8-diazaspiro[4.5]decan-1-one`), `1` unknown (`5.0%`), `2` invalid (`10.0%`), `0` prices (PubChem publishes vendor IDs without prices).

---

## 3. Model

`buy-or-make` combines deterministic cheminformatics normalization, tiered stock/retrosynthesis search, conservative unit-aware cost algebra, and optional pretrained synthetic-complexity scoring:

1. **Structure Normalization & Desalting (`buyormake/canonical.py`)**:
   - Parses SMILES via `Chem.MolFromSmiles`, strips reaction atom-map tags, caches `num_heavy_atoms` and parsed `Chem.Mol`, and optionally desalts/uncharges counterions (`rdMolStandardize.FragmentParent` + `Uncharger`) so pharmaceutical HCl/sodium salts match neutral catalog stocks transparently.
2. **Stock Lookup Backends (`buyormake/stock.py`)**:
   - `LocalStockBackend`: Ingests `.tsv`, `.csv`, `.jsonl`, `.smi`/`.smi.gz`, priced/unpriced `.sqlite` (`PRAGMA table_info(stock)`), and `.hdf5` (`pandas`/PyTables or `h5py` fallback, raising `RuntimeError` if unreadable).
   - `PubChemVendorBackend`: Resolves `InChIKey -> CID` (`GET` + throttled SMILES `POST` fallback), fetches `Chemical Vendors` from PUG-View (`max_offers=25`), retries `429`/`5xx`/`TransportError` with exponential backoff, never caches transient failures as negative results, and supports configurable cache TTL (`default=30` days online, infinite offline).
3. **Multi-Tier & `Retro*` Best-First Retrosynthesis Engine (`buyormake/routes.py`, `buyormake/retro_star.py`)**:
   - **Frontier Planner — `RetroStarPlanner` (`planner="retro_star"`)**: Implements Chen et al. (ICML 2020) best-first AND-OR tree search up to `max_steps=6`, propagating target-rooted synthesis cost numbers $V_t(m \mid \mathcal{T})$ and admissible minimum reaction-step bounds across multi-branch convergent and linear trees, guided by calibrated `SCScore` + `SAscore` value heuristics $V_0(m)$, cycle-free ancestor InChIKey tracking, `100` atom-conserved medicinal chemistry SMARTS templates (`EXTENDED_RETRO_STAR_TEMPLATES`), and exact mixed-SKU / shared-synthon discrete pack-size MOQ outlay minimization (`optimize_discrete_pack_order`, `optimize_route_moq_outlay`).
   - **Tier 1 — `AiZynthFinderBackend`** (optional `[aizynthfinder]` extra): Runs MCTS up to `max_steps`, preserves MCTS `Route.score`, selects solved branches across multi-OR nodes, and falls through in `auto` mode if no route has all-purchasable leaves.
   - **Tier 2 — `ReactionLookupBackend`**: Multi-step lookup against locally indexed `PaRoutes` (`CC-BY-4.0`), `USPTO` (`CC0-1.0` / `MIT`), `ORDerly` (`CC-BY-4.0`), and `ORD` (`CC-BY-SA-4.0`) reactions in `data/reactions.sqlite`, with `parse_reaction_smiles_for_indexing` separating catalysts, solvents, and inorganic spectators into `agents_json`.
   - **Tier 3 — `TemplateFallbackBackend` (`planner="lookup_smarts"`)**: Built-in 1- and 2-step RDKit `ReactionFromSmarts` engine over `17` core templates (`15` acyclic `!@` bimolecular disconnections + `2` monocyclic `@` lactone/lactam ring-openings), enforcing exact per-element heavy-atom conservation against each template's leaving-group delta (`RETRO_TEMPLATE_DELTAS`) and chemoselectivity guards.
4. **Stock Availability, Unit-Safe Cost Bases & Strict-Null Budget Gate (`buyormake/cost.py`, `buyormake/rank.py`, `buyormake/handoff.py`)**:
   - **Stock-Flag Honesty (`purchasable_offers`)**: Filters out any `Offer` explicitly marked `in_stock=False` across candidate ranking (`evaluate_candidates`), leaf purchasability (`_resolve_leaf_offers`, `_is_purchasable`), route cost/MOQ optimization, and `labshift` scheduling handoff (`to_procurement_constraint` / `build_labshift_release_map`), while preserving `in_stock=None` for free registries. When every offer for a molecule is `in_stock=False`, `RankResult` preserves `best_price` as evidence, sets `status="unknown"` and `all_offers_out_of_stock=True`, and `to_procurement_constraint` sets `schedulable=False` and `release_time_minutes=None`.
   - **Distinct Cost Field Bases (`<currency>/g` vs. `<currency>`)**:
     - `best_price` (when mass `pack_size` is known), `est_cost`, `lower_bound_cost`, and `est_cost_yield_adjusted` (`RouteCostBreakdown.yield_adjusted_cost`) are intensive **per-1-gram** quantities labelled `<currency>/g` (e.g. `"USD/g"`), invariant to `target_scale_g`:
       $$\text{est\_cost}\;[\text{USD/g}] = \sum_{\ell \in \text{unique leaves}} \left(\min_{o \in \text{purchasable}(\ell)} \text{price/g}(o)\right) \cdot \sum_{d \in \text{depths}(\ell)} \frac{\text{stoich\_excess}}{\text{step\_yield}^d} \cdot \frac{\text{MW}_\ell}{\text{MW}_{\text{target}}}$$
     - `min_outlay_cost` is an extensive **whole-order discrete pack cash outlay** for `target_scale_g` grams labelled in bare `<currency>` (e.g. `"USD"`): $\sum_\ell \min_{o} \lceil m_{\ell}(\text{target\_scale\_g}) / \text{pack\_g}(o) \rceil \cdot \text{price}(o)$.
     - Both `buyable` (unless `allow_raw_pack_budget=True`) and `makeable` rows require `<budget_currency>/g` units under a per-gram `budget`, rejecting incompatible units with `price_status="incompatible_unit"`. When `enforce_moq_outlay=True` and `max_outlay_budget` is omitted, the discrete pack cash-outlay limit scales as `budget * target_scale_g`.
   - **Strict-Null Budget Default (`require_complete_cost`)**: Whenever `budget` or `max_outlay_budget` is set, `require_complete_cost` defaults to `True` in `filter_decisions_from_ranked`, `filter_proposals`, and CLI `filter`, rejecting unpriced (`best_price=None`) and incomplete-cost (`cost_complete=False`) rows by default in exact agreement with `BuyOrMakeLean.unpriced_never_passes_budget` (explicit opt-out via `require_complete_cost=False` / `--allow-incomplete-cost`).
5. **Lean 4.36.0 Formal Guardrails & Differential Parity (`lean/BuyOrMakeLean.lean`, `scripts/check_lean_axioms_and_parity.py`)**:
   - Verified in CI with `0 sorry`, `0` custom axioms, and **175 `#eval` differential parity cases** (`85` directly against shipping `buyormake.*` code + `90` against `community_extensions.pareto_green_triage`):
     - `checkBuyOrMake`, `checkBuyOrMake_ok_iff`, `unpriced_never_passes_budget`, `non_usd_per_gram_rejected_under_budget`: checked directly against shipping `buyormake.rank.evaluate_candidates` + `buyormake.rank.filter_decisions_from_ranked` (`shipping_gate`, 40 local-catalog cases) and `community_extensions.pareto_green_triage.verify_option_guardrails` (`proposals`, 60 cases, delegating budget checks to `filter_decisions_from_ranked`).
     - `routeStoichiometryOK`, `ring_cleavage_rejected`, `atom_imbalance_rejected`, `cross_halogen_swap_rejected`, `suzuki_boron_imbalance_rejected`, `main_group_or_metal_imbalance_rejected`: checked against `buyormake.routes.TemplateFallbackBackend.one_step`, `buyormake.retro_star.RetroStarPlanner.expand_extended_smarts`, and `atom_vec_from_smiles`.
     - `discretePackCount`, `discrete_pack_covers_required_mass`, `shared_synthon_pack_count_subadditive`, `shared_synthon_outlay_subadditive`, `moqOutlayPolicyOK`, `moq_outlay_gate_sound`: checked directly against shipping `buyormake.cost.min_outlay_for_offers` (`moq`, 45 cases).
     - `compoundYieldMassUg`, `yield_compounded_ge_lower_bound`, `yield_compounded_monotone_depth`, `paretoDominates`, `paretoDominates_irrefl`: checked against `community_extensions.pareto_green_triage.pareto_dominates` and fixed-point step-depth yield compounding (`yield_pareto`, 30 cases).
6. **Synthetic Accessibility Scorers (`buyormake/scores.py`)**:
   - `SAscore` (`Ertl & Schuffenhauer`, `rdkit.Contrib.SA_Score`, range `1`–`10`).
   - `SCScore` (`Coley et al.`, standalone NumPy MLP over `1024`-bit chiral radius-2 Morgan fingerprints mapped via overflow-safe $1 + 4 \cdot \tfrac{1}{2}(1 + \tanh(z/2)) \in [1, 5]$, returning `1.0` floor on zero-heavy-atom inputs).

---

## 4. Training

- **Zero Training Performed in Repository**: `buy-or-make` does **not** train or fine-tune any neural network weights.
- **Frozen Pretrained / Heuristic Components**:
  - `SAscore` uses RDKit's deterministic fragment-frequency lookup table plus spiro/bridgehead/macrocycle/chiral complexity penalties.
  - `SCScoreModel` (`buyormake/scores.py`) caches and runs pure NumPy forward inference over Coley et al.'s MIT-licensed Reaxys-pretrained checkpoint (`model.ckpt-10654.as_numpy.json.gz`, SHA-256 `fa90bd552307b2efd05eb1e53033bc9eadadfd5bdca682e541040982af99536b`).
  - Optional `AiZynthFinder` uses AstraZeneca's CC-BY-4.0 USPTO-pretrained ONNX expansion/filter policy weights (`scripts/fetch_aizynth.py`).

---

## 5. History

| Commit / Branch | Date | Author | Summary |
|---|---|---|---|
| `76991c8` (`main`) | `2026-11-18` | Yani Guan | Initial release (`Add buy-or-make: cost-aware availability check for candidate molecules`): `buyormake` package, 5 exported JSON Schemas, 8 fetch scripts, 20-molecule CPU demo, and `pytest` suite. |
| `bbfa649` (`origin/review-fixes`) | `2026-11-18` | Yani Guan | Set `cost_complete = False` on `buyable` rows, prevented `cheapest_offer_for_target` from comparing mixed currencies or mismatched non-mass pack sizes, added `estimate_route_cost_detailed` / `route_cost_method`, marked ZINC indirect license check, strengthened `SCScore` test, and clarified demo zero-price caveats. |
| `fix/resolve-audit-issues` (`cc2d33e`) | `2026-11-22` | Yani Guan | Enforced acyclic `!@` SMARTS + exact element conservation (`Analog-6` true SnAr disconnection), guarded PubChem cache against transient `5xx`/`429`/timeout errors with throttled `POST` + retries + TTL, propagated route/leaf licenses, partitioned `_rank_sort_key` and `filter_proposals` unit/completeness policies, added multi-step lookup + spectator filtering + AiZynth score/OR-node/stock fixes, added salt fallback + InChIKey deduplication, priced SQLite + `h5py` HDF5 support, atomic SHA-256 downloads, CLI `filter`/`--records-out`/`.smi` names, and `assayfold` + `labshift` handoff helpers. |
| `fix/phase2-audit-and-wire` (`#32`, `5fd1afb`) | `2026-11-25` | Yani Guan | Resolved Phase 2 audit issues (`#26`–`#31`: symmetric desalting, optional stereo-relaxed 14-char skeletal disclosure, MOQ outlay & yield-adjusted cost bounds, chemoselectivity SMARTS guards, SQLite WAL/readonly locking, and `BuyOrMakeLean` proofs) and wired `autolab_wire.py`. |
| `feat/open-datasets-2026` (`#33`, `3be1f4a`) | `2026-11-25` | Yani Guan | Audited and integrated 7 open 2024–2026 chemical stock, compound, reaction, and retrosynthesis corpora (`PaRoutes 2.0`, `ChEMBL 37`, `SureChEMBL`, `COCONUT 2.0`, `ORDerly`, `Syntheseus USPTO-50K & Retro*`, `SynPlanner`) with `Source.spdx_license` schema fields, `buyormake/open_corpora.py` checksummed multi-format adapters, and `scripts/fetch_open_datasets_2026.py`. |
| `feat/frontier-retrostar-paroutes-and-lean-ci` (`#34`, `ec33997`) | `2026-11-27` | Yani Guan | Added Chen et al. (ICML 2020) `Retro*` best-first AND-OR search (`buyormake/retro_star.py`), explicit `data_source` provenance separating `synthetic_reference_catalog` from vendor/stopping stocks, fail-closed `autolab_wire.py` validation, Lean 4.36.0 CI, and `240`-target `PaRoutes 2.0`/`SynPlanner` + `500`-reaction held-out `USPTO-50K` benchmark. |
| `fix/audit-35-38-leakfree-stock-units-lean` | `2026-11-28` | Yani Guan | Resolved audit Issues `#35`–`#38`: eliminated target reference-route leakage in `scripts/benchmark_retrosynthesis.py` (whole-tree + product InChIKey/connectivity-block exclusion, `index_hygiene`, `100` atom-conserved `Retro*` SMARTS templates, all `19/19` corpus SHA-256/MD5 digests pinned), enforced `Offer.in_stock=False` via `purchasable_offers` across ranking/routes/`labshift` handoff (`all_offers_out_of_stock`), labelled route per-gram costs `USD/g` vs whole-order MOQ `USD` (`budget * target_scale_g`), defaulted `require_complete_cost=True` under active budgets, and added the `40`-case `shipping_gate` Lean differential parity suite (`175` total cases). |

---

## 6. Unique Properties

1. **Strict Epistemic Honesty on Prices & Reachability**: Never imputes synthetic catalog prices (`price` stays `null`, `price_basis` stays `"unknown"` when unlisted) and never conflates "not found within `max_steps`" (`"unknown"`) with synthetic impossibility.
2. **Atom-Conserved Retrosynthesis & Stoichiometric Lower Bounds**: Enforces acyclic `!@` bond constraints plus exact heavy-atom delta checks on every SMARTS disconnection and strips unmapped patent/ORD spectators before stoichiometric cost estimation.
3. **End-to-End License Gate Propagation**: Aggregates target, reaction-database (`CC-BY-SA-4.0` ORD, `CC0-1.0` USPTO, `CC-BY-4.0` AiZynth), leaf-stock, and scorer licenses into `AvailabilityRecord.licenses` and selects the most restrictive governing license when sealing `AssayfoldRecord`s.
4. **Seamless `L5` -> `L0`/`L5` Auto-Lab Integration**: Converts triaged molecules directly into content-addressed [`assayfold`](https://github.com/yaniguan/assayfold) records (`to_assayfold_record()`) and [`labshift`](https://github.com/yaniguan/labshift) job release-time constraints (`to_procurement_constraint()` / `build_labshift_release_map()`).

---

## 7. 15-Repository Auto-Lab Wire Contract (`autolab-wire/v1`)

`autolab_wire.py` (`REPO_ID = "buy-or-make"`, `LAYER = "L5"`) implements the standardized `autolab-wire/v1` packet interface connecting `buy-or-make` across the 15-repository Auto-Lab closed loop:

- **Upstream (`L6` / `L0`)**: Accepts `candidate_smiles` and `budget_usd_per_g` after epistemic hypothesis screening ([`lab-ask-why`](https://github.com/yaniguan/lab-ask-why)) and split-conformal trust-region gating ([`lab-decide`](https://github.com/yaniguan/lab-decide)), while maintaining lossless molecule, offer, and provenance field compatibility with [`assayfold`](https://github.com/yaniguan/assayfold) (`AvailabilityRecord.to_assayfold_record()`).
- **Downstream (`L5` / `L4`)**: Emits `out["procurement"]` (`smiles`, `status`, `listed_price_usd_per_g`, `n_steps`, `price_basis`, `within_budget`, `lead_time_days`, `lead_time_complete`, `release_time_minutes`, `moq_pack_g`, `selected_pack_size`, `min_order_outlay_usd`, `precursor_smiles`, `precursor_basis`, `reaction_steps`) to [`labshift`](https://github.com/yaniguan/labshift) and [`lab-attrib`](https://github.com/yaniguan/lab-attrib) prior to digital-twin rehearsal (`lab-rehearsal`) and capability-gated execution (`lab-verify` -> `stepbot`).
- **Offline Determinism & Lean Guardrails**: Resolves repository paths exclusively via `Path(__file__).resolve().parent` (never `os.getcwd()`), executes strictly offline against a deterministic built-in `listed_catalog` reference stock (aspirin `"CC(=O)Oc1ccccc1C(=O)O"` at `$0.42/g` plus standard synthons), preserves `None` prices (`status="unknown"`, `within_budget=False`) for unlisted molecules, and enforces `BuyOrMakeLean`'s `unpriced_never_passes_budget` and acyclic `!@` stoichiometric conservation invariants (`halt_reason="L5_PROCUREMENT_BUDGET_BLOCK"` when `require_priced_within_budget=True`).
- **Full 15-Repository Ecosystem**: `lab-ears` (`L2`), `lab-ask-why` (`L6`), `lab-decide` (`L6`), `lab-automation-atlas` (`L5`), `buy-or-make` (`L5`), `labshift` (`L5`), `lab-rehearsal` (`L4`), `lab-verify` (`L4`), `stepbot` (`L1`), `lab-hands` (`L1`), `lab-feet` (`L0`), `lab-attrib` (`L4` post-execution causal blame, NNLS precursor deconvolution & SHA-256 certified ledger feeding closed-loop recalibration back into `lab-decide` and `labbench-posttrain`), `twinstep` (`L3`), `labfill` (`L3`), and `labbench-posttrain` (`L3`).
