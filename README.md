# buy-or-make

Check whether candidate molecules are directly purchasable from chemical vendors or makeable in N steps from purchasable building blocks, with conservative leaf-only cost lower bounds computed strictly from real listed catalog prices. Missing prices stay `null` and unverified molecules stay `"unknown"`.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"

# Optional: fetch MIT-licensed SCScore NumPy weights (6.18 MB, SHA-256 verified) into data/scscore/
python scripts/fetch_scscore.py

# Rank or budget-filter candidates from a SMILES file (SMILES [optional_name])
cat > /tmp/candidates.smi << 'EOF'
CC(=O)Oc1ccccc1C(=O)O Aspirin
COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1 Analog-3
invalid_smiles BadMol
EOF

buy-or-make rank /tmp/candidates.smi --max-steps 2
buy-or-make rank /tmp/candidates.smi --max-steps 2 --records-out /tmp/records.jsonl --json
buy-or-make filter /tmp/candidates.smi --max-steps 2 --budget 25.0 --budget-currency USD --require-complete-cost --json
```

> **Auto-Lab Ecosystem**: Layer `L5 · Executive Planning, Procurement & Scheduling` in [auto-lab](https://github.com/yaniguan/auto-lab) · [Project Studio](https://yaniguan.github.io/auto-lab/)

## System Profile: Object · Data · Model · Training · History · Unique Properties

Full architectural documentation, cost/ranking invariants, commit and branch evolution, and sibling-repository contracts are maintained in [`docs/ARCHITECTURE_AND_HISTORY.md`](docs/ARCHITECTURE_AND_HISTORY.md):

- **Object**: Cost-aware molecular availability and retrosynthetic procurement triage (`L5`) verifying whether candidate molecules from `L6`/`L5` planners (`lab-ask-why`, `lab-decide`) are directly purchasable or makeable in $N$ steps from purchasable building blocks before wet-lab scheduling (`labshift`).
- **Data**: 19 license-audited chemical vendor, building-block, natural-product, patent-compound, reaction, and retrosynthesis benchmark sources (`buyormake/sources.py` + `buyormake/open_corpora.py`, gitignored `data/` fetch-only isolation with atomic `.part` checksum verification), 5 strict Pydantic v2 / JSON Schema contracts (`extra="forbid"`, aligned with `assayfold`), and a 20-molecule CPU benchmark (`examples/demo.py`).
- **Model**: RDKit canonical/InChIKey normalization with optional salt/counterion parent fallback, live throttled/retry-hardened `PubChem PUG-REST/PUG-View` + multi-format `LocalStockBackend` (including priced SQLite and `h5py`/PyTables HDF5), 3-tier retrosynthesis (`AiZynthFinder` MCTS, multi-step `USPTO`/`ORD` SQLite lookup with spectator filtering, and 17 atom-conserved RDKit SMARTS templates), unit-aware `USD/g` stoichiometric leaf cost estimator, and `SAscore` + overflow-safe NumPy `SCScore`.
- **Training**: Zero model training in-repo; combines deterministic cheminformatics rules with RDKit's fragment `SAscore` and frozen pretrained Reaxys `SCScore` weights (`6.18 MB`, optional `CC-BY-4.0` `AiZynthFinder` ONNX policies).
- **History**: Initialized on `main` in `76991c8` (`2026-11-18`), refined in `origin/review-fixes` (`bbfa649`, `2026-11-18`) for `cost_complete`, single-currency/same-pack price comparability, excluded-leaf reporting, and indirect ZINC license-check markers, and hardened on `fix/resolve-audit-issues` with acyclic `!@` SMARTS + element conservation, transient HTTP cache guards, route/leaf license propagation, unit-aware ranking/budget gates, and `assayfold`/`labshift` handoff adapters.
- **Unique Properties**: Strict `null`/`"unknown"` honesty (never invents missing catalog prices or conflates search bounds with synthetic impossibility), completeness-aware leaf-cost lower bounds, zero-redistribution license gating across target/route/leaf sources, and deterministic SHA-256 `AvailabilityRecord` + `AssayfoldRecord` sealing compatible with `assayfold` and `labshift`.

## Repository Layout

```text
buy-or-make/
├── .github/workflows/
│   └── ci.yml                            # Matrix pytest + ruff + zero-schema-drift CI (Python 3.10–3.13)
├── buyormake/
│   ├── __init__.py                       # Public API exports
│   ├── cache.py                          # Persistent SQLite API response cache with TTL & refresh
│   ├── canonical.py                      # Canonical JSON, SHA-256, atomic downloads & RDKit desalting
│   ├── cli.py                            # CLI entry point (rank, filter, export-schema, sources)
│   ├── cost.py                           # Mass-normalized USD/g offer & leaf-only route cost logic
│   ├── handoff.py                        # L5 -> L0/L5 adapters (assayfold Record & labshift release times)
│   ├── open_corpora.py                   # Checksummed 2024-2026 open stock, compound & reaction adapters
│   ├── rank.py                           # cost_aware_rank, evaluate_candidates, filter_proposals
│   ├── routes.py                         # AiZynthFinder, multi-step USPTO/ORD lookup & 17 SMARTS templates
│   ├── schema.py                         # Pydantic v2 models (extra='forbid') & JSON Schema exporter
│   ├── scores.py                         # RDKit SAscore & overflow-safe standalone NumPy SCScore
│   ├── sources.py                        # 19 license-audited data & model Source definitions
│   └── stock.py                          # Multi-format LocalStockBackend & PubChemVendorBackend
├── community_extensions/                 # 3D cost-lead-PMI Pareto router & shared-synthon cart optimizer
├── docs/
│   └── ARCHITECTURE_AND_HISTORY.md       # Deep technical profile (Object, Data, Model, History, Wire)
├── examples/
│   ├── demo.py                           # 20-molecule cold vs. cached CPU benchmark script
│   └── demo_output.md                    # Recorded CPU runtime and candidate triage table
├── lean/                                 # Lean 4 formalization (BuyOrMakeLean guardrails & theorems)
├── schema/                               # Exported JSON Schemas (offer, route, source, rank_result, ...)
├── scripts/                              # Checksum-verified local fetchers & reaction/SDF indexers
├── tests/                                # Offline pytest suite (wire, cache, cli, cost, handoff, ...)
├── autolab_wire.py                       # Layer L5 Auto-Lab 15-repo wire node (autolab-wire/v1)
├── LICENSE                               # MIT License
├── pyproject.toml                        # Build metadata, dependencies, and optional extras
└── README.md                             # Usage guide, source license table, demo metrics & caveats
```

## Python API

```python
from buyormake import (
    build_labshift_release_map,
    cost_aware_rank,
    evaluate_candidates,
    filter_proposals,
)

candidates = [
    "CC(=O)Oc1ccccc1C(=O)O",                    # buyable
    "COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1",    # novel analog
    "C1CC(invalid",                              # invalid
]

ranked, records = evaluate_candidates(
    candidates,
    max_steps=2,
    backends=("local_stock", "pubchem"),
    route_backend="auto",
    offline=False,
)
for r in ranked:
    print(r.status, r.n_steps, r.best_price, r.est_cost, r.reasons[0])

# Convert to sealed assayfold Record(kind='molecule') or labshift job release offsets
af_rec = records[0].to_assayfold_record()
release_minutes = build_labshift_release_map(records)

# Co-scientist proposal filter (strict-null USD/g budget gate by default; pass require_complete_cost=False to opt out)
summary = filter_proposals(
    candidates,
    max_steps=2,
    budget=25.0,
    budget_currency="USD",
)
print("kept:", len(summary.kept), "over_budget:", len(summary.over_budget))
```

Ranking order: `buyable` (`USD/g` ascending, then other units grouped by unit, then `null` price) < `makeable` (fewer steps, complete `USD/g` cost before other currencies before partial cost before `null`) < `unknown` < `invalid`.

## 15-Repository Auto-Lab Wire Contract (`autolab-wire/v1`)

[`autolab_wire.py`](autolab_wire.py) exposes the deterministic offline `L5` wire node (`REPO_ID = "buy-or-make"`, `LAYER = "L5"`) for the 15-repository Auto-Lab loop:

- **Upstream & Downstream Links**: Consumes candidate molecules admitted by [`lab-decide`](https://github.com/yaniguan/lab-decide) (`L6` calibrated conformal decision gate) with molecule/offer schemas aligned to [`assayfold`](https://github.com/yaniguan/assayfold) (`L0` open assay hub), and hands off verified procurement status and lead-time constraints to [`labshift`](https://github.com/yaniguan/labshift) (`L5` multi-instrument window scheduler).
- **Packet Semantics (`process_packet(packet: dict) -> dict`)**:
  - Resolves all repository paths via `Path(__file__).resolve().parent` (never `os.getcwd()`) and runs strictly offline (zero live HTTP requests) against a deterministic built-in `listed_catalog` reference stock (e.g., aspirin `"CC(=O)Oc1ccccc1C(=O)O"` at `$0.42/g` and standard building blocks).
  - Writes `out["procurement"]` containing core triage fields (`smiles`, `status`, `listed_price_usd_per_g`, `n_steps`, `price_basis`, `within_budget`) alongside `labshift` scheduling/MOQ fields (`lead_time_days`, `lead_time_complete`, `release_time_minutes`, `moq_pack_g`, `selected_pack_size`, `min_order_outlay_usd`) and `lab-attrib` stoichiometric attribution records (`precursor_smiles`, `precursor_basis`, `reaction_steps`), keeping unlisted prices strictly `None` (`price_basis="unknown"`, `within_budget=False`).
  - Enforces the Lean 4 `unpriced_never_passes_budget` and acyclic `!@` heavy-atom conservation invariants (`lean/BuyOrMakeLean.lean`): when `require_priced_within_budget=True` and `listed_price_usd_per_g is None` or exceeds `budget_usd_per_g`, sets `halted=True` and `halt_reason="L5_PROCUREMENT_BUDGET_BLOCK"`.
  - Appends structured execution records to `out["trace"]` and `out["bottom_reasoning_traces"]` (`invariant="unpriced_never_passes_budget"`).
- **15-Repository Stack**: Operates alongside `lab-ears` (`L2`), `lab-ask-why` (`L6`), `lab-decide` (`L6`), `lab-automation-atlas` (`L5`), `buy-or-make` (`L5`), `labshift` (`L5`), `lab-rehearsal` (`L4`), `lab-verify` (`L4`), `stepbot` (`L1`), `lab-hands` (`L1`), `lab-feet` (`L0`), `lab-attrib` (`L4` post-execution causal blame, NNLS XRD deconvolution & SHA-256 certified ledger), `twinstep` (`L3`), `labfill` (`L3`), and `labbench-posttrain` (`L3`).

## Schema summary

All models in [`buyormake/schema.py`](buyormake/schema.py) use Pydantic v2 (`extra="forbid"`), enforce finite non-negative numeric invariants, align shared field names (`Source`, `Quantity`, `Uncertainty`, `Provenance`, `Molecule`, `origin`, `schema_version`, `record_hash`) with the `assayfold` hub schema, and export JSON Schemas to `schema/`:

- `Source`: evaluated dataset/API metadata, `license` + explicit `spdx_license` (SPDX ID or `LicenseRef-Custom-Terms`), `license_evidence`, `license_quote`, ISO `license_checked` (`YYYY-MM-DD`), `license_check_method` (`"direct" | "indirect"`), `commercial_use_ok`, `redistribution_ok`, `status`.
- `Offer`: `inchikey`, `vendor`, `catalog_id`, `pack_size` (`Quantity`), `purity`, `price` (`float | null`), `currency` (`str | null`), `price_basis` (`"listed_catalog" | "unknown"`), `lead_time_days`, `in_stock`, `source_id`, `source_url`, `fetched_utc`, `license`.
- `Route` / `RouteStep`: `target_inchikey`, `n_steps`, `tree_depth`, `steps`, `leaves`, `leaves_purchasable`, `method`, `score`, `est_cost` (`Quantity | null`, per `1 g` target in `<currency>/g`), `est_cost_yield_adjusted` (`Quantity | null`, per `1 g` target in `<currency>/g`), `min_outlay_cost` (`Quantity | null`, discrete pack MOQ cash outlay for the full `target_scale_g` order in `<currency>`), `est_cost_method`, `est_cost_complete`.
- `AvailabilityRecord`: `molecule`, `offers`, `routes`, `sa_scores` (`sascore`, `scscore`, `rascore`), `origin`, `provenance` (propagating target, route-step, and leaf sources), `licenses`, `schema_version`, deterministic `record_hash`, and `.to_assayfold_record()`.
- `RankResult`: `candidate`, `name`, `canonical_smiles`, `inchikey`, `status` (`buyable | makeable | unknown | invalid`), `n_steps`, `best_price`, `est_cost`, `est_cost_yield_adjusted`, `min_outlay_cost`, `target_scale_g`, `cost_complete`, `all_offers_out_of_stock`, `sa_scores`, `reasons`, `evidence`.

Regenerate JSON Schemas with:

```bash
buy-or-make export-schema --out schema
```

## Evaluated sources & licenses

Every source license/terms page was checked directly (`2026-11-18` for the initial 12 sources and `2026-11-25` for the 7 audited 2024–2026 open chemistry corpora), except ZINC, whose terms page was unreachable and was checked indirectly (`license_checked = "2026-11-18"`, `license_check_method = "indirect"` in `buyormake/sources.py`). Non-redistributable datasets and large bulk dumps are never committed to the repository; fetch scripts stream via `.part` temporary files with checksum verification into the gitignored `data/` directory.

| Source | Role | SPDX / License | Commercial use | Redistribution | Status | Evidence |
|---|---|---|---|---|---|---|
| PubChem PUG-REST & PUG-View | Live chemical vendors & catalog IDs (no prices) | `LicenseRef-Custom-Terms` (`custom/terms`) | Yes (per depositor) | No (cached locally only) | `integrated` | [NCBI Policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) |
| ZINC (ZINC15/ZINC20/ZINC22) | Purchasable / in-stock SMILES subsets | `LicenseRef-Custom-Terms` (`custom/terms`) | Yes ("free to use for everyone") | No (no major portions without written permission) | `fetch_only` | [ZINC Terms](https://wiki.docking.org/index.php/Terms_And_Conditions) |
| Enamine Building Blocks | Building block SDFs (manual registration) | `LicenseRef-Custom-Terms` (`custom/terms`) | Yes for own use after registration ("Open Structure Data" section) | No (no public republishing) | `fetch_only` | [Enamine Terms](https://enamine.net/terms-of-use) |
| eMolecules Free Download | Free SMILES + IDs (350 MB; prices require paid ePlus) | `LicenseRef-Custom-Terms` (`custom/terms`) | No (full commercial requires ePlus) | No | `fetch_only` | [eMolecules Readme](https://downloads.emolecules.com/free/2026-11-01/Readme.txt) |
| Mcule Database | Free in-stock & BB SMILES + IDs (no prices) | `LicenseRef-Custom-Terms` (`custom/terms`) | Yes (own use) | No | `fetch_only` | [Mcule Terms](https://mcule.com/terms/) |
| AiZynthFinder Code & Models | MCTS retrosynthesis + USPTO policy models | `CC-BY-4.0` (models) / `MIT` (code) | Yes | Yes | `fetch_only` | [AiZynthFinder LICENSE](https://raw.githubusercontent.com/MolecularAI/aizynthfinder/master/LICENSE) |
| AiZynthFinder ZINC Stock | `zinc_stock_17_04_20.hdf5` (663 MB, ZINC-derived) | `LicenseRef-Custom-Terms` (`custom/terms`) | No | No | `fetch_only` | [Figshare 12334577](https://figshare.com/articles/dataset/12334577) |
| USPTO Reactions (Lowe) | Reported patent reactions (1976–Sep 2016) | `CC0-1.0` | Yes | Yes | `fetch_only` | [Figshare 5104873](https://api.figshare.com/v2/articles/5104873) |
| Open Reaction Database (`ord-data`) | Reported reactions (parquet shards on Hugging Face) | `CC-BY-SA-4.0` | Yes | Yes (ShareAlike) | `fetch_only` | [ord-data LICENSE](https://raw.githubusercontent.com/open-reaction-database/ord-data/main/LICENSE) |
| RDKit Contrib SA_Score | Synthetic accessibility heuristic (1–10) | `BSD-3-Clause` | Yes | Yes | `integrated` | [RDKit license.txt](https://raw.githubusercontent.com/rdkit/rdkit/master/license.txt) |
| SCScore (Coley et al.) | Synthetic complexity score (1–5, NumPy) | `MIT` | Yes | Yes | `integrated` | [SCScore LICENSE](https://raw.githubusercontent.com/connorcoley/scscore/master/LICENSE) |
| RAscore (Reymond Group) | Learned retrosynthetic accessibility (52 MB) | `MIT` | Yes | Yes | `not_used` | [RAscore LICENSE](https://raw.githubusercontent.com/reymond-group/RAscore/master/LICENSE) |
| PaRoutes 2.0 (`paroutes`) | Multi-step retrosynthesis routes (`n1`/`n5`), stocks & reactions | `CC-BY-4.0` (data) / `Apache-2.0` (code) | Yes | Yes | `fetch_only` | [Zenodo 7341155](https://zenodo.org/api/records/7341155) |
| ChEMBL 37 (`chembl`) | Bioactive compound SMILES & InChIKeys (May 2026 release) | `CC-BY-SA-3.0` | Yes | Yes (ShareAlike) | `fetch_only` | [ChEMBL 37 README](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/README) |
| SureChEMBL (`surechembl`) | Patent chemistry Parquet corpus (`bulk_data/latest`, Nov 2026) | `CC-BY-4.0` | Yes | Yes | `fetch_only` | [SureChEMBL Docs](https://chembl.gitbook.io/surechembl/llms-full.txt) |
| COCONUT 2.0 (`coconut`) | ~700k open natural product structures (`2026-11` dump) | `CC0-1.0` | Yes | Yes | `fetch_only` | [COCONUT About](https://coconut.naturalproducts.net/about) |
| ORDerly (`orderly`) | Cleaned ORD forward/retrosynthesis Parquet benchmarks (JCIM 2024) | `CC-BY-4.0` | Yes | Yes | `fetch_only` | [Figshare 23298467](https://api.figshare.com/v2/articles/23298467) |
| Syntheseus USPTO-50K & Retro* (`syntheseus_uspto_50k`) | Single-step `USPTO-50K` & `Retro*` stock inventory (2024–2026) | `MIT` | Yes | Yes | `fetch_only` | [Figshare 25459573](https://api.figshare.com/v2/articles/25459573) |
| SynPlanner Open Data (`synplanner_data`) | Standardized USPTO reactions, building blocks & SAscore targets | `MIT` | Yes | Yes | `fetch_only` | [HF SynPlanner](https://huggingface.co/api/models/Laboratoire-De-Chemoinformatique/SynPlanner) |

Note on ZINC check: `wiki.docking.org` dropped TLS connections from the test host on `2026-11-18`; its terms (`fetch_only`, no major-portion redistribution without written permission) were verified from prior documentation and search results. This is machine-visible on the `Source` record via `license_check_method = "indirect"` and `notes` starting with `LICENSE CHECK INDIRECT`.

## Multi-Step Retrosynthesis & Single-Step Benchmark (`PaRoutes 2.0`, `SynPlanner`, `USPTO-50K`)

<!-- BENCHMARK_TABLE_START -->
Measured offline on CPU (`Python 3.11.15`, `RDKit 2026.03.4`, seed `20261127`) via [`scripts/benchmark_retrosynthesis.py`](scripts/benchmark_retrosynthesis.py) -> [`results/retrosynthesis_benchmark.json`](results/retrosynthesis_benchmark.json) over **240 stratified multi-step targets** (`PaRoutes 2.0 n1` `n=80`, `PaRoutes 2.0 n5` `n=80`, `SynPlanner` SAscore-stratified `n=80`) plus **500 held-out `Syntheseus USPTO-50K` test reactions** (`raw_test.csv` strictly excluded). **Split & target-exclusion hygiene**: the headline leakage-free reaction index (`106,599` reactions across `19,840` non-target `PaRoutes` trees + `USPTO-50K` train+val) excludes all `160` sampled `PaRoutes` target trees **and** any reaction whose product `InChIKey` or 14-char connectivity block matches any of the `240` evaluated targets (`all_target_product_inchikey_overlap = 0`):

#### Overall & Dataset-Stratified Multi-Step Solve Rates (`n = 240`)

| Planner / Configuration | Target Routes Indexed? | Overall Solve % | `PaRoutes n1` | `PaRoutes n5` | `SynPlanner` | Leaf Stock Completeness | Median Search Time |
|---|:---:|---:|---:|---:|---:|---:|---:|
| 1. Stock-only lookup (`max_steps=0`) | **No (0 overlap)** | 0.0% (0/240) | 0.0% | 0.0% | 0.0% | 0.0% (solved: 100.0%) | 1.07 ms |
| 2. 1-step held-out lookup + 2-step greedy SMARTS (`17` templates, leakage-free) | **No (0 overlap)** | 2.1% (5/240) | 2.5% | 0.0% | 3.8% | 34.4% (solved: 100.0%) | 18.08 ms |
| 3a. Leakage-free `Retro*` (`max_steps=4`, held-out `PaRoutes` + `USPTO-50K` + `100` SMARTS) | **No (0 overlap)** | 27.9% (67/240) | 31.2% | 41.2% | 11.2% | 49.0% (solved: 100.0%) | 294.35 ms |
| 3b. **Headline — Leakage-free `Retro*`** (`max_steps=5`, held-out `PaRoutes` + `USPTO-50K` + `100` SMARTS) | **No (0 overlap)** | 30.0% (72/240) | 35.0% | 43.8% | 11.2% | 50.4% (solved: 100.0%) | 322.37 ms |
| 3c. Strict zero-`PaRoutes` `Retro*` (`max_steps=5`, `USPTO-50K` train+val + `100` SMARTS only) | **No (0 overlap)** | 24.6% (59/240) | 26.2% | 36.2% | 11.2% | 47.3% (solved: 100.0%) | 326.20 ms |
| 4. *Reference-route retrieval upper bound* (`max_steps=5`, target's own `PaRoutes` route indexed — NOT headline) | Yes (upper bound) | 70.4% (169/240) | 100.0% | 100.0% | 11.2% | 78.8% (solved: 100.0%) | 340.62 ms |

#### Solve Rate Stratified by `PaRoutes 2.0` Reference Route Depth & `SAscore` Band

| Planner | Depth 2 | Depth 3 | Depth 4 | Depth 5 | Convergent Subset | `SA [1.5,2.5)` | `SA [2.5,3.5)` | `SA [3.5,4.5)` | `SA [4.5,10]` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1. Stock-only lookup (`max_steps=0`) | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| 2. 1-step held-out lookup + 2-step greedy SMARTS (`17` templates, leakage-free) | 0.0% | 5.0% | 0.0% | 0.0% | 0.0% | 4.6% | 0.9% | 0.0% | 4.5% |
| 3a. Leakage-free `Retro*` (`max_steps=4`, held-out `PaRoutes` + `USPTO-50K` + `100` SMARTS) | 45.0% | 45.0% | 27.5% | 27.5% | 40.0% | 46.2% | 26.1% | 13.2% | 9.1% |
| 3b. **Headline — Leakage-free `Retro*`** (`max_steps=5`, held-out `PaRoutes` + `USPTO-50K` + `100` SMARTS) | 47.5% | 47.5% | 27.5% | 35.0% | 50.0% | 49.2% | 28.7% | 13.2% | 9.1% |
| 3c. Strict zero-`PaRoutes` `Retro*` (`max_steps=5`, `USPTO-50K` train+val + `100` SMARTS only) | 47.5% | 32.5% | 20.0% | 25.0% | 30.0% | 41.5% | 23.5% | 7.9% | 9.1% |
| 4. *Reference-route retrieval upper bound* (`max_steps=5`, target's own `PaRoutes` route indexed — NOT headline) | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 80.0% | 83.5% | 47.4% | 13.6% |

#### Held-Out `Syntheseus USPTO-50K` Single-Step Exact Reactant-Set Accuracy (`n = 500`)

| Single-Step Disconnection Engine | Top-1 | Top-3 | Top-5 | Top-10 | Median Time |
|---|---:|---:|---:|---:|---:|
| `USPTO-50K` train+val lookup + `17` core SMARTS | 11.4% | 16.0% | 16.6% | 16.6% | 5.61 ms |
| `USPTO-50K` train+val lookup + `100` `Retro*` SMARTS | 17.6% | 36.6% | 44.6% | 49.4% | 25.75 ms |

#### Verified Input Dataset SHA-256 Digests

| Artifact | License | Bytes | SHA-256 Digest |
|---|---|---:|---|
| `stock_n1.txt` (`paroutes_stock_n1`) | `CC-BY-4.0` | 368,802 | `b8641e2028846d995953b30f7acd36e7521d0301159518d172497340a0bfb115` |
| `stock_n5.txt` (`paroutes_stock_n5`) | `CC-BY-4.0` | 360,970 | `199ec2d9fad0638c3ae53f2ea091dcd6bf9207a938b12ae1aa9429c37dd097af` |
| `ref_routes_n1.json` (`paroutes_routes_n1`) | `CC-BY-4.0` | 64,187,987 | `b818f2115ce10da51f84bf13f19dc2206863b1b074ef348017de258e9e1bff02` |
| `ref_routes_n5.json` (`paroutes_routes_n5`) | `CC-BY-4.0` | 82,120,928 | `d878bec737d34a861b859912ac783f065c380f704f6d0087af4ea8ef667fa2bf` |
| `uspto50k.zip` (`syntheseus_uspto_50k`) | `MIT` | 19,940,280 | `0f3871191a47f234ba5f5f48fca7fea08169718fe6d057c998a7588d8372f5c8` |
| `building_blocks_em_sa_ln.smi` (`synplanner_building_blocks`) | `MIT` | 5,784,208 | `53b6d5eec7b3c522e8499213b95f5884028ce73fd83e09e1e63b6ded97e5affb` |
| `targets_with_sascore_1.5_2.5.smi` (`synplanner_targets_sa_1_5_2_5`) | `MIT` | 4,361 | `b3e7263bd803a2dbe694d0a8a8d7562971ab9fc6a1994d306104eec5fe7eca28` |
| `targets_with_sascore_2.5_3.5.smi` (`synplanner_targets_sa_2_5_3_5`) | `MIT` | 4,875 | `68bc68a98549646e5509c2e46135f318517e4e96905575c9e38e30303160a61e` |
| `targets_with_sascore_3.5_4.5.smi` (`synplanner_targets_sa_3_5_4_5`) | `MIT` | 5,398 | `d8125674587b9b2ecfcc268ace17adc8f1cf35c2281a001cbc9487d744d42856` |
| `targets_with_sascore_4.5_5.5.smi` (`synplanner_targets_sa_4_5_5_5`) | `MIT` | 6,087 | `9b02a61ee21876b1977ac319d4305e9cfb201be531041953caceda517182516b` |
| `model.ckpt-10654.as_numpy.json.gz` (`scscore_weights`) | `MIT` | 6,178,760 | `fa90bd552307b2efd05eb1e53033bc9eadadfd5bdca682e541040982af99536b` |
<!-- BENCHMARK_TABLE_END -->

## Demo coverage and CPU runtime

Measured from an actual run of [`examples/demo.py`](examples/demo.py) on `2026-11-22` (20 candidates: 6 common drugs, 6 reagents/building blocks, 6 hand-drawn analogs, 2 invalid SMILES; 2 of the 6 analogs are listed directly by PubChem vendors, and 3 are 1-step makeable via atom-conserved acyclic SMARTS disconnections) using live PubChem PUG-REST/PUG-View + RDKit SAscore + NumPy SCScore + built-in 2-step template fallback:

- **Hardware / Environment**: AuthenticAMD CPU, Python 3.11.15, RDKit 2026.03.4
- **Cold run wall time**: `38.737 s` total (`1.937 s / molecule`, throttled on both `GET` and `POST` at <= 5 req/s to PubChem)
- **Cached run wall time**: `0.058 s` total (`0.003 s / molecule`, offline SQLite cache)

| Metric | Count | Fraction of total (n=20) | Fraction of valid (n=18) |
|---|---:|---:|---:|
| Buyable (listed by >=1 vendor in PubChem) | 14 | 70.0% | 77.8% |
| Makeable (complete route <= 2 steps) | 3 | 15.0% | 16.7% |
| Unknown (not in stock; no complete route <= 2 steps) | 1 | 5.0% | 5.6% |
| Invalid SMILES | 2 | 10.0% | - |
| Complete price known | 0 | 0.0% | 0.0% |
| Route found (for non-buyable) | 3 | 15.0% | 16.7% |

The demo configuration has no price source: PubChem Chemical Vendors publishes vendors and catalog IDs but not prices, and no local priced offer file was supplied. The cost columns (`best_price`, `est_cost`) are therefore structurally `null`, and the "Complete price known = 0" row is a consequence of that configuration, not a measured result of the cost logic. The cost estimator (`cheapest_offer_for_target`, `estimate_route_cost`) is demonstrated by the unit tests in `tests/test_cost.py`, `tests/test_rank.py`, and `tests/test_stock_formats.py`, which use synthetic local offers.

Full candidate-by-candidate table is saved in [`examples/demo_output.md`](examples/demo_output.md).

## Optional backends & local data scripts

- **Local priced/unpriced stock** (`.tsv`, `.csv`, `.jsonl`, `.smi`, `.smi.gz`, `.sqlite`, `.hdf5`): pass `--stock path/to/file` or place `data/stock.sqlite` / `data/stock.tsv`. Install `pip install -e ".[hdf5]"` for HDF5 stocks.
- **USPTO / ORD exact reaction lookup** (supports up to 3-step chaining and filters unmapped catalyst/solvent spectators):
  ```bash
  pip install -e ".[data]"
  python scripts/fetch_uspto.py --max-rows 100000
  python scripts/fetch_ord.py --max-files 1
  ```
- **AiZynthFinder (optional heavy backend, ~790 MB models/stock + 178 packages)**:
  ```bash
  pip install -e ".[aizynthfinder]"
  python scripts/fetch_aizynth.py --with-stock
  buy-or-make rank candidates.smi --route-backend aizynthfinder
  ```
- **2024–2026 Open Corpora (`PaRoutes`, `ChEMBL 37`, `SureChEMBL`, `COCONUT 2.0`, `ORDerly`, `Syntheseus USPTO-50K`, `SynPlanner`)**:
  ```bash
  python scripts/fetch_open_datasets_2026.py --list
  python scripts/fetch_open_datasets_2026.py --artifact paroutes_stock_n1
  python scripts/fetch_open_datasets_2026.py --artifact syntheseus_uspto_50k --max-rows 50000
  python scripts/fetch_open_datasets_2026.py --artifact chembl_37_chemreps --max-rows 100000
  ```
  Stopping-stock artifacts import into `data/stock.sqlite`, reaction/route benchmarks index into `data/reactions.sqlite`, and non-vendor reference corpora (`ChEMBL`, `SureChEMBL`, `COCONUT 2.0`) import into `data/reference_compounds.sqlite` (pass `--stock data/reference_compounds.sqlite` only when deliberately matching against reference compounds).
- **ZINC / eMolecules / Mcule / Enamine**: see [`scripts/fetch_zinc.py`](scripts/fetch_zinc.py), [`scripts/fetch_emolecules.py`](scripts/fetch_emolecules.py), [`scripts/fetch_mcule.py`](scripts/fetch_mcule.py), and [`scripts/fetch_enamine.py`](scripts/fetch_enamine.py).

## Caveats and limitations

- **Not purchasing advice.** Catalog prices change frequently and vary by region, currency, purity, pack size, institutional discounts, and account terms. A publicly listed price is not a binding quote.
- **PubChem vendor lists are capped at 25 per molecule** in the default backend (`max_offers=25`). Online cache entries default to a 30-day TTL (`--cache-max-age-days`, bypassable via `--refresh-cache`), while `--offline` uses cached entries without expiration. Transient HTTP `5xx`/`429`/timeout errors are retried and never cached as negative results.
- **Missing prices stay `null`.** Free public sources (PubChem Chemical Vendors, eMolecules free, Mcule free, ZINC, Enamine Open Structure Data) do not publish prices. Unless you supply a local offer file with real `listed_catalog` prices, `best_price` and `est_cost` remain `null`.
- **Explicitly out-of-stock listings do not make a molecule buyable or schedulable.** Offers with `in_stock=False` (`purchasable_offers(offers) == []`) never qualify a candidate as `buyable` or a retrosynthetic leaf as purchasable (`leaves_purchasable=False`). When a molecule has catalog offers that are all `in_stock=False`, `evaluate_candidates` retains the offers and `best_price` as evidence, sets `status="unknown"` and `all_offers_out_of_stock=True` with a reason naming the out-of-stock vendors, and `to_procurement_constraint` sets `schedulable=False` and `release_time_minutes=None` (omitting the molecule from `build_labshift_release_map`). Offers with `in_stock=None` (free registries that publish no stock flag) remain purchasable as before.
- **Cost field bases (`<currency>/g` per-gram vs. `<currency>` whole-order cash outlay).** On `Route` and `RankResult`, `best_price` (when mass `pack_size` is known), `est_cost`, `lower_bound_cost`, and `est_cost_yield_adjusted` are intensive **per-1-gram** quantities labelled `<currency>/g` (e.g. `USD/g`, invariant to `target_scale_g`), whereas `min_outlay_cost` is an extensive **whole-order discrete pack cash outlay** for `target_scale_g` grams labelled `<currency>` (e.g. `USD`). When `enforce_moq_outlay=True` and no explicit `max_outlay_budget` is passed, `filter_decisions_from_ranked` derives the cash-outlay limit as `budget * target_scale_g` rather than reusing the per-gram `budget` directly.
- **Route cost is a rough lower bound.** When leaf prices are supplied, `est_cost` sums `cheapest_USD_per_g * count * (MW_leaf / MW_target)` per 1 g target (`unit="USD/g"`). It ignores solvents, reagents, catalysts, reaction yield, purification, analytical QC, and labor. If any leaf lacks a price, `cost_complete` is `false`. Leaves left out of the sum (unparseable SMILES, out of stock, or no listed USD/g price) are listed in the route's `est_cost_method` and in the `reasons` of a makeable result. For buyable molecules `est_cost` is `null` and `cost_complete` is `false`; the listed price is in `best_price`.
- **`est_cost` for makeable molecules is the cost of a found route, not the minimum-cost route.** The template search is not exhaustive: it carries at most 6 first-step disconnections into a second step, tries at most 6 second-step disconnections for each, and stops expanding a first step once one second step reaches purchasable leaves. The reported route is the best of those enumerated (fewest steps, then cost); a cheaper route may exist.
- **Prices in different currencies or pack sizes are never compared, and active budgets reject unpriced candidates by default.** `best_price` is the cheapest per-gram price in `USD/g` when mass pack sizes are known; otherwise raw pack prices are compared only when all offers share one pack size, and non-USD offers only when they all use a single currency. In cross-molecule ranking and `filter_proposals(budget=...)`, `USD/g` offers and `USD/g` route costs are never numerically interleaved with raw per-pack `USD` or non-USD offers (`price_status="incompatible_unit"`). Whenever a numeric `budget` (or `max_outlay_budget`) is active, `require_complete_cost` defaults to `True` (matching `BuyOrMakeLean.unpriced_never_passes_budget`; pass `require_complete_cost=False` or `--allow-incomplete-cost` to opt out).
- **Template fallback is heuristic.** The 17 core RDKit SMARTS disconnections (`routes.py`, max 2 steps) and 100 extended medicinal-chemistry `Retro*` SMARTS templates (`retro_star.py`, verified for exact heavy-atom conservation and chemoselectivity guards) do not replace experimental validation of reaction conditions or stereocontrol.
- **SA/SC scores are heuristics**, not synthesis guarantees.
- **`"unknown"` does not mean impossible.** It only means no stock offer or complete route within `max_steps` was found in the configured backends.

## License & author

- Code: MIT License ([LICENSE](LICENSE)).
- Data & model weights: each source retains its own license (see table above).
- Author: Yani Guan.
