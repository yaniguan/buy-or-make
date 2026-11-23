# Architecture and History — `buy-or-make`

> **Auto-Lab Ecosystem Role**: Layer `L5 · Executive Planning, Procurement & Scheduling` in [Robot AI Auto-Lab (`auto-lab`)](https://github.com/yaniguan/auto-lab) ([Project Studio](https://yaniguan.github.io/auto-lab/)).

---

## 1. Object

### Problem Statement
Generative chemistry models, active-learning loops, and LLM co-scientists frequently propose candidate molecules without checking whether they can actually be bought off the shelf or synthesized in a small number of steps from commercially available building blocks—or worse, rely on hallucinated prices and opaque synthetic accessibility heuristics as proxies for real procurement cost.

`buy-or-make` solves **cost-aware molecular procurement and retrosynthetic triage**: given a list of candidate SMILES strings, it normalizes chemical structures (`RDKit`), queries local stock catalogs and live vendor registries (`PubChem PUG-REST / PUG-View`), searches for short retrosynthetic disconnections (`1`–`N` steps) back to purchasable leaves when a molecule is not directly in stock, computes conservative leaf-only stoichiometric cost lower bounds strictly from real `listed_catalog` prices, and ranks or filters proposals under explicit budget constraints while keeping missing prices `null` and unverified compounds `"unknown"`.

### Role in the Auto-Lab `L0`–`L6` Hierarchy
Within the [Auto-Lab](https://github.com/yaniguan/auto-lab) Brain/Body architecture (`L0`–`L6`), `buy-or-make` operates in **`L5 · Executive Planning, Procurement & Scheduling`**:
- **Upstream (`L6` / `L5` / `L3` Brain)**: Filters candidate molecule batches proposed by hypothesis generators ([`lab-ask-why`](https://github.com/yaniguan/lab-ask-why)), trust-region Bayesian optimizers ([`lab-decide`](https://github.com/yaniguan/lab-decide)), or molecular extraction/design workflows ([`ouroboros-ocsr`](https://github.com/yaniguan/ouroboros-ocsr), [`YLM`](https://github.com/yaniguan/YLM), [`WignerFlow`](https://github.com/yaniguan/WignerFlow)).
- **Peer (`L5`)**: Works alongside [`labshift`](https://github.com/yaniguan/labshift) (multi-instrument scheduling once reagents are sourced) and [`lab-automation-atlas`](https://github.com/yaniguan/lab-automation-atlas) (ecosystem map of automation platforms and datasets).
- **Downstream (`L0` Ledger / `L4` Verification)**: Exports Pydantic v2 `AvailabilityRecord` objects sealed with deterministic SHA-256 `record_hash` digests using field names (`Source`, `Quantity`, `Uncertainty`, `AgentRef`, `SoftwareRef`, `Provenance`, `Molecule`, `origin`, `schema_version`, `record_hash`) directly aligned with the [`assayfold`](https://github.com/yaniguan/assayfold) (`L0`) open-data hub schema.

### Inputs, Outputs, and Scope Boundaries
- **Inputs**:
  - Candidate SMILES strings (via Python API `cost_aware_rank` / `evaluate_candidates` / `filter_proposals` or CLI `buy-or-make rank <file.smi>`).
  - Optional local stock files (`.tsv`, `.csv`, `.jsonl`, `.smi`, `.smi.gz`, `.sqlite`, `.hdf5`), reaction lookup SQLite index (`data/reactions.sqlite`), SCScore weights (`data/scscore/`), AiZynthFinder config (`data/aizynth/config.yml`), step cap `max_steps`, and optional `budget`.
- **Outputs**:
  - Sorted `list[RankResult]` classifying each candidate into `buyable | makeable | unknown | invalid`, with `n_steps`, `best_price` (`Quantity | null`), `est_cost` (`Quantity | null`), `cost_complete`, `sa_scores` (`SAScores`), human-readable `reasons`, and vendor `evidence` URLs.
  - Input-ordered `list[AvailabilityRecord]` sealed with canonical SHA-256 `record_hash`.
  - `FilterSummary` (`kept`, `over_budget`, `rejected`, per-candidate `ProposalDecision`) for co-scientist budget gating.
  - Five strict JSON Schema files (`additionalProperties: false`) under `schema/`.
- **Scope Boundaries**:
  - **Not purchasing quotes**: Catalog prices vary by pack size, purity, region, currency, and institutional contracts; free public sources (PubChem, ZINC, eMolecules free, Mcule free, Enamine Open Structure Data) publish zero prices, so `best_price` and `est_cost` remain `null` unless local `listed_catalog` prices are supplied.
  - **Lower-bound stoichiometry only**: Route `est_cost` sums $1\text{ equiv}$ leaf building-block mass per $1\text{ g}$ target and excludes solvents, catalysts, reagents, non-unity yields, purification, QC, and labor.

---

## 2. Data

### Evaluated Source Registry & License Gate (`buyormake/sources.py`)
Twelve external chemical vendor, reaction, stock, and scoring sources were audited on `2026-11-18` and encoded as validated `Source` objects (`buy-or-make sources`):

| Source ID | Modality / Role | License | Commercial OK | Redistribution OK | Status |
|---|---|---|---|---|---|
| `pubchem_pug` | Live vendor names, RegistryIDs & URLs (`<= 5 req/s`; no prices) | `custom/terms` (NCBI + depositor terms) | Yes | **No** (local SQLite cache only) | `integrated` |
| `rdkit_sascore` | Fragment + complexity heuristic (`1`–`10`) | `BSD-3-Clause` | Yes | Yes | `integrated` |
| `scscore` | Pretrained neural synthetic complexity (`1`–`5`, `6.18 MB`) | `MIT` | Yes | Yes | `integrated` |
| `zinc` | Purchasable SMILES tranches (`scripts/fetch_zinc.py`) | `custom/terms` | Yes | **No** | `fetch_only` |
| `enamine_bb` | Building block SDFs after registration (`scripts/fetch_enamine.py`) | `custom/terms` | Yes (own use) | **No** | `fetch_only` |
| `emolecules` | Free `version.smi.gz` (`~350 MB`, no prices; `scripts/fetch_emolecules.py`) | `custom/terms` | **No** (requires ePlus) | **No** | `fetch_only` |
| `mcule` | Free in-stock building blocks (`~151 MB`; `scripts/fetch_mcule.py`) | `custom/terms` | Yes (own use) | **No** | `fetch_only` |
| `aizynthfinder_models` | MCTS retrosynthesis + USPTO ONNX models (`scripts/fetch_aizynth.py`) | `CC-BY-4.0` (models) / `MIT` (code) | Yes | Yes | `fetch_only` |
| `aizynthfinder_zinc_stock` | `zinc_stock_17_04_20.hdf5` (`663 MB`, ZINC-derived) | `custom/terms` | **No** | **No** | `fetch_only` |
| `uspto_lowe` | Lowe 1976–Sep 2016 patent reactions (`75.2 MB`; `scripts/fetch_uspto.py`) | `CC0-1.0` | Yes | Yes | `fetch_only` |
| `ord_data` | Open Reaction Database Parquet shards (`scripts/fetch_ord.py`) | `CC-BY-SA-4.0` | Yes | Yes (ShareAlike) | `fetch_only` |
| `rascore` | Reymond Group retrosynthetic classifier (`52.4 MB`) | `MIT` | Yes | Yes | `not_used` |

All fetched files, SQLite indices (`data/stock.sqlite`, `data/reactions.sqlite`), weights (`data/scscore/`), and HTTP caches (`data/cache.sqlite`) live inside gitignored `data/` paths (`.gitignore`) so non-redistributable vendor records are never committed.

### Unified Pydantic v2 & JSON Schemas (`buyormake/schema.py`, `schema/`)
All domain models inherit `_Model(BaseModel)` with `ConfigDict(extra="forbid")`, rejecting undeclared or hallucinated fields at validation time. Canonical hashing (`AvailabilityRecord.content_hash()`) serializes sorted-key, whitespace-free UTF-8 JSON (`allow_nan=False`), strips `record_hash` and volatile `provenance.recorded_utc`, and computes a 64-char SHA-256 digest.

### 20-Molecule Benchmark Demo (`examples/demo.py`, `examples/demo_output.md`)
Evaluates 20 curated inputs (`6` common drugs, `6` building blocks/reagents, `6` novel analogs, `2` malformed SMILES) on CPU (`Intel Xeon`, Python `3.12.16`, RDKit `2026.09.2`):
- **Cold live PubChem run**: `24.970 s` (`1.249 s / molecule` under `0.22 s` request spacing).
- **Warm offline SQLite cache run**: `0.099 s` (`0.005 s / molecule`, `~252x` speedup).
- **Triage breakdown**: `14` buyable (`70.0%`), `3` makeable (`15.0%`), `1` unknown (`5.0%`), `2` invalid (`10.0%`), `0` prices (since PubChem publishes vendor listings without prices).

---

## 3. Model

`buy-or-make` combines deterministic cheminformatics normalization, tiered stock/retrosynthesis search, conservative unit-aware cost algebra, and optional pretrained synthetic-complexity scoring:

1. **Structure Normalization (`buyormake/canonical.py`)**:
   - Parses SMILES via `Chem.MolFromSmiles`, generates canonical isomeric SMILES (`MolToSmiles(..., isomericSmiles=True, canonical=True)`), standard InChIKey (`Chem.MolToInchiKey`), molecular formula, and molecular weight (`Descriptors.MolWt`).
2. **Stock Lookup Backends (`buyormake/stock.py`)**:
   - `LocalStockBackend`: Ingests `.tsv`, `.csv`, `.jsonl`, `.smi`/`.smi.gz`, `.sqlite`, and AiZynthFinder `.hdf5` stocks keyed by InChIKey.
   - `PubChemVendorBackend`: Resolves `InChIKey -> CID` via PUG-REST (with SMILES `POST` fallback), fetches `Chemical Vendors` depositions from PUG-View (`max_offers=25`), throttles live calls (`min_interval_sec=0.22`), and persists responses in `ResponseCache` (`data/cache.sqlite`).
3. **Three-Tier Retrosynthesis Engine (`buyormake/routes.py`)**:
   - **Tier 1 — `AiZynthFinderBackend`** (optional `[aizynthfinder]` extra): Runs Monte Carlo Tree Search (MCTS) using USPTO ONNX expansion/filter policy networks up to `max_transforms = max_steps`.
   - **Tier 2 — `ReactionLookupBackend`**: Exact `product_inchikey` lookup against locally indexed `USPTO` (`CC0-1.0`) and `ORD` (`CC-BY-SA-4.0`) reactions in `data/reactions.sqlite`.
   - **Tier 3 — `TemplateFallbackBackend`**: Built-in 1- and 2-step heuristic RDKit `ReactionFromSmarts` engine over `14` curated reaction classes (`amide_coupling_secondary/tertiary`, `esterification`, `sulfonamide_secondary/tertiary`, `urea_formation`, `carbamate_formation`, `suzuki_biaryl`, `suzuki_aryl_alkenyl`, `williamson_aryl_ether`, `williamson_alkyl_ether`, `reductive_amination_secondary/tertiary`, `buchwald_hartwig_cn`), capped at depth `2` and `max_branch_per_node=12`, preferring balanced fragment sizes.
4. **Conservative Cost Estimator & Lexicographic Ranker (`buyormake/cost.py`, `buyormake/rank.py`)**:
   - Normalizes mass-based offers (`ug`, `µg`, `mg`, `g`, `kg`) to `USD/g`.
   - Estimates 1-equiv stoichiometric leaf cost per $1\text{ g}$ target:
     $$\text{est\_cost} = \sum_{\ell \in \text{unique leaves}} \left(\min_{o \in \text{offers}(\ell)} \text{USD/g}(o)\right) \cdot \text{count}(\ell) \cdot \frac{\text{MW}_\ell}{\text{MW}_{\text{target}}}$$
   - Marks `est_cost_complete = True` iff *every* unique leaf has a valid `USD/g` price.
   - Ranks candidates strictly by `buyable` (known price ascending, then `null` price) $<$ `makeable` (`n_steps` ascending, complete cost before partial before `null`, cost ascending) $<$ `unknown` $<$ `invalid`.
5. **Synthetic Accessibility Scorers (`buyormake/scores.py`)**:
   - `SAscore` (`Ertl & Schuffenhauer`, bundled with `rdkit.Contrib.SA_Score`, range `1`–`10`).
   - `SCScore` (`Coley et al.`, standalone pure-NumPy 5-layer ReLU MLP over `1024`-bit chiral radius-2 Morgan fingerprints mapped via $1 + 4\,\sigma(z) \in [1, 5]$).

---

## 4. Training

- **Zero Training Performed in Repository**: `buy-or-make` does **not** train or fine-tune any neural network weights.
- **Frozen Pretrained / Heuristic Components**:
  - `SAscore` uses RDKit's deterministic fragment-frequency lookup table plus spiro/bridgehead/macrocycle/chiral complexity penalties.
  - `SCScoreModel` (`buyormake/scores.py`) performs pure NumPy forward inference over Coley et al.'s MIT-licensed Reaxys-pretrained checkpoint (`model.ckpt-10654.as_numpy.json.gz`, fetched on demand via `scripts/fetch_scscore.py`; gracefully returns `scscore = None` if absent).
  - Optional `AiZynthFinder` uses AstraZeneca's CC-BY-4.0 USPTO-pretrained ONNX expansion/filter policy weights (`scripts/fetch_aizynth.py`).

---

## 5. History

| Commit | Date | Author | Branch / Ref | Summary |
|---|---|---|---|---|
| `76991c8` (`76991c8b668a28a7e82454ebf8c122cf2e4f73e3`) | `2026-11-18` | Yani Guan | `main` | Initial release (`Add buy-or-make: cost-aware availability check for candidate molecules`): implemented `buyormake` package (`canonical`, `schema`, `sources`, `cache`, `stock`, `cost`, `routes`, `scores`, `rank`, `cli`), 5 exported JSON Schemas (`schema/`), 8 data/model fetch scripts (`scripts/`), 20-molecule CPU demo (`examples/`), 7 `pytest` modules (`tests/`), and Python 3.10–3.12 CI (`.github/workflows/ci.yml`). |
| `bbfa649` (`bbfa649148c0410d414735e49e03a3e7f34d0eac`) | `2026-11-18` | Yani Guan | `origin/review-fixes` | Open review-fixes branch (`Fix cost_complete, mixed-currency prices, and demo claims`): sets `cost_complete = False` on `buyable` rows (reserving `cost_complete` strictly for route cost estimates), prevents `cheapest_offer_for_target` from comparing mixed currencies or mismatched non-mass pack sizes, adds `estimate_route_cost_detailed` / `route_cost_method` to surface excluded leaves in `est_cost_method` and `RankResult.reasons`, marks ZINC's indirect license verification machine-visibly (`CHECKED_DATE_INDIRECT`), strengthens the 2-layer hand-verified `SCScore` unit test, and clarifies demo zero-price configuration caveats in `README.md`. |

---

## 6. Unique Properties

1. **Strict Epistemic Honesty on Prices & Reachability**: Never imputes synthetic catalog prices (`price` stays `null`, `price_basis` stays `"unknown"` when unlisted) and never conflates "not found within `max_steps`" (`"unknown"`) with synthetic impossibility.
2. **Conservative Leaf-Only Lower Bounds with Completeness Flags**: Separates complete stoichiometric leaf sums (`est_cost_complete = True`) from partial lower bounds (`est_cost_complete = False`), ranking complete estimates ahead of incomplete ones at equal step count.
3. **Zero-Redistribution License-Gated Data Architecture**: Encodes explicit `commercial_use_ok`, `redistribution_ok`, and `license_quote` provenance on all 12 evaluated sources while isolating downloaded vendor catalogs and API caches inside gitignored `data/`.
4. **Deterministic Content-Addressed `assayfold` Compatibility**: Every `AvailabilityRecord` shares core value types (`Source`, `Quantity`, `Uncertainty`, `Provenance`, `Molecule`, `record_hash`) with [`assayfold`](https://github.com/yaniguan/assayfold) (`L0`) and seals records with timestamp-invariant canonical SHA-256 hashes, enabling lossless ingestion into the Auto-Lab experiment ledger and downstream scheduling in [`labshift`](https://github.com/yaniguan/labshift).
