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

### Evaluated Source Registry & License Gate (`buyormake/sources.py`)
Twelve external chemical vendor, reaction, stock, and scoring sources were audited on `2026-11-18` and encoded as validated `Source` objects (`buy-or-make sources`, with ISO `YYYY-MM-DD` `license_checked` and explicit `license_check_method="direct"|"indirect"`):

| Source ID | Modality / Role | License | Commercial OK | Redistribution OK | Status |
|---|---|---|---|---|---|
| `pubchem_pug` | Live vendor names, RegistryIDs & URLs (`<= 5 req/s`; no prices) | `custom/terms` (NCBI + depositor terms) | Yes | **No** (local SQLite cache only) | `integrated` |
| `rdkit_sascore` | Fragment + complexity heuristic (`1`–`10`) | `BSD-3-Clause` | Yes | Yes | `integrated` |
| `scscore` | Pretrained neural synthetic complexity (`1`–`5`, `6.18 MB`) | `MIT` | Yes | Yes | `integrated` |
| `zinc` | Purchasable SMILES tranches (`scripts/fetch_zinc.py`, `indirect` check) | `custom/terms` | Yes | **No** | `fetch_only` |
| `enamine_bb` | Building block SDFs after registration (`scripts/fetch_enamine.py`) | `custom/terms` | Yes (own use) | **No** | `fetch_only` |
| `emolecules` | Free `version.smi.gz` (`~350 MB`, no prices; `scripts/fetch_emolecules.py`) | `custom/terms` | **No** (requires ePlus) | **No** | `fetch_only` |
| `mcule` | Free in-stock building blocks (`~151 MB`; `scripts/fetch_mcule.py`) | `custom/terms` | Yes (own use) | **No** | `fetch_only` |
| `aizynthfinder_models` | MCTS retrosynthesis + USPTO ONNX models (`scripts/fetch_aizynth.py`) | `CC-BY-4.0` (models) / `MIT` (code) | Yes | Yes | `fetch_only` |
| `aizynthfinder_zinc_stock` | `zinc_stock_17_04_20.hdf5` (`663 MB`, ZINC-derived) | `custom/terms` | **No** | **No** | `fetch_only` |
| `uspto_lowe` | Lowe 1976–Sep 2016 patent reactions (`75.2 MB`; `scripts/fetch_uspto.py`) | `CC0-1.0` | Yes | Yes | `fetch_only` |
| `ord_data` | Open Reaction Database Parquet shards (`scripts/fetch_ord.py`) | `CC-BY-SA-4.0` | Yes | Yes (ShareAlike) | `fetch_only` |
| `rascore` | Reymond Group retrosynthetic classifier (`52.4 MB`) | `MIT` | Yes | Yes | `not_used` |

All fetched files stream into `<dest>.part`, verify pinned SHA-256 / MD5 / LFS digests via `buyormake.canonical.atomic_download`, and atomically `os.replace` into gitignored `data/` paths (`.gitignore`).

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
3. **Three-Tier Retrosynthesis Engine (`buyormake/routes.py`)**:
   - **Tier 1 — `AiZynthFinderBackend`** (optional `[aizynthfinder]` extra): Runs MCTS up to `max_steps`, preserves MCTS `Route.score`, selects solved branches across multi-OR nodes, and falls through in `auto` mode if no route has all-purchasable leaves.
   - **Tier 2 — `ReactionLookupBackend`**: Bounded `1`–`3` step lookup against locally indexed `USPTO` (`CC0-1.0`) and `ORD` (`CC-BY-SA-4.0`) reactions in `data/reactions.sqlite`, with `parse_reaction_smiles_for_indexing` separating catalysts, solvents, and inorganic spectators into `agents_json`.
   - **Tier 3 — `TemplateFallbackBackend`**: Built-in 1- and 2-step RDKit `ReactionFromSmarts` engine over `17` templates (`15` acyclic `!@` bimolecular disconnections + `2` monocyclic `@` lactone/lactam ring-openings), enforcing exact per-element heavy-atom conservation against each template's leaving-group delta (`RETRO_TEMPLATE_DELTAS`).
4. **Unit-Safe Cost Estimator, Ranker & Budget Gate (`buyormake/cost.py`, `buyormake/rank.py`)**:
   - Normalizes mass-based offers (`ug`, `µg`, `mg`, `g`, `kg`) to `<currency>/g`.
   - Estimates 1-equiv stoichiometric leaf cost per $1\text{ g}$ target:
     $$\text{est\_cost} = \sum_{\ell \in \text{unique leaves}} \left(\min_{o \in \text{offers}(\ell)} \text{price/g}(o)\right) \cdot \text{count}(\ell) \cdot \frac{\text{MW}_\ell}{\text{MW}_{\text{target}}}$$
   - Separates `buyable` sort tiers into normalized `USD/g`, other currencies/raw-pack units grouped by `unit`, and `null` price; enforces unit compatibility and `require_complete_cost` in `filter_proposals`.
5. **Synthetic Accessibility Scorers (`buyormake/scores.py`)**:
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
| `fix/resolve-audit-issues` | `2026-11-22` | Yani Guan | Enforced acyclic `!@` SMARTS + exact element conservation (`Analog-6` true SnAr disconnection), guarded PubChem cache against transient `5xx`/`429`/timeout errors with throttled `POST` + retries + TTL, propagated route/leaf licenses, partitioned `_rank_sort_key` and `filter_proposals` unit/completeness policies, added multi-step lookup + spectator filtering + AiZynth score/OR-node/stock fixes, added salt fallback + InChIKey deduplication, priced SQLite + `h5py` HDF5 support, atomic SHA-256 downloads, CLI `filter`/`--records-out`/`.smi` names, and `assayfold` + `labshift` handoff helpers. |

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
