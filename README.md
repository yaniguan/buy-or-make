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
- **Data**: 12 license-audited chemical vendor, reaction, stock, and scoring sources (`buyormake/sources.py`, gitignored `data/` fetch-only isolation with atomic `.part` checksum verification), 5 strict Pydantic v2 / JSON Schema contracts (`extra="forbid"`, aligned with `assayfold`), and a 20-molecule CPU benchmark (`examples/demo.py`).
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
│   ├── rank.py                           # cost_aware_rank, evaluate_candidates, filter_proposals
│   ├── routes.py                         # AiZynthFinder, multi-step USPTO/ORD lookup & 17 SMARTS templates
│   ├── schema.py                         # Pydantic v2 models (extra='forbid') & JSON Schema exporter
│   ├── scores.py                         # RDKit SAscore & overflow-safe standalone NumPy SCScore
│   ├── sources.py                        # 12 license-audited data & model Source definitions
│   └── stock.py                          # Multi-format LocalStockBackend & PubChemVendorBackend
├── docs/
│   └── ARCHITECTURE_AND_HISTORY.md       # Deep technical profile (Object, Data, Model, History)
├── examples/
│   ├── demo.py                           # 20-molecule cold vs. cached CPU benchmark script
│   └── demo_output.md                    # Recorded CPU runtime and candidate triage table
├── schema/                               # Exported JSON Schemas (offer, route, source, rank_result, ...)
├── scripts/                              # Checksum-verified local fetchers & reaction/SDF indexers
├── tests/                                # Offline pytest suite (cache, cli, cost, handoff, indexers, ...)
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

# Co-scientist proposal filter (unit-checked per-gram budget + optional completeness enforcement)
summary = filter_proposals(
    candidates,
    max_steps=2,
    budget=25.0,
    budget_currency="USD",
    require_complete_cost=False,
)
print("kept:", len(summary.kept), "over_budget:", len(summary.over_budget))
```

Ranking order: `buyable` (`USD/g` ascending, then other units grouped by unit, then `null` price) < `makeable` (fewer steps, complete `USD` cost before other currencies before partial cost before `null`) < `unknown` < `invalid`.

## Schema summary

All models in [`buyormake/schema.py`](buyormake/schema.py) use Pydantic v2 (`extra="forbid"`), enforce finite non-negative numeric invariants, align shared field names (`Source`, `Quantity`, `Uncertainty`, `Provenance`, `Molecule`, `origin`, `schema_version`, `record_hash`) with the `assayfold` hub schema, and export JSON Schemas to `schema/`:

- `Source`: evaluated dataset/API metadata, SPDX or `custom/terms`, `license_evidence`, `license_quote`, ISO `license_checked` (`YYYY-MM-DD`), `license_check_method` (`"direct" | "indirect"`), `commercial_use_ok`, `redistribution_ok`, `status`.
- `Offer`: `inchikey`, `vendor`, `catalog_id`, `pack_size` (`Quantity`), `purity`, `price` (`float | null`), `currency` (`str | null`), `price_basis` (`"listed_catalog" | "unknown"`), `lead_time_days`, `in_stock`, `source_id`, `source_url`, `fetched_utc`, `license`.
- `Route` / `RouteStep`: `target_inchikey`, `n_steps`, `steps`, `leaves`, `leaves_purchasable`, `method`, `score`, `est_cost` (`Quantity | null`), `est_cost_method`, `est_cost_complete`.
- `AvailabilityRecord`: `molecule`, `offers`, `routes`, `sa_scores` (`sascore`, `scscore`, `rascore`), `origin`, `provenance` (propagating target, route-step, and leaf sources), `licenses`, `schema_version`, deterministic `record_hash`, and `.to_assayfold_record()`.
- `RankResult`: `candidate`, `name`, `canonical_smiles`, `inchikey`, `status` (`buyable | makeable | unknown | invalid`), `n_steps`, `best_price`, `est_cost`, `cost_complete`, `sa_scores`, `reasons`, `evidence`.

Regenerate JSON Schemas with:

```bash
buy-or-make export-schema --out schema
```

## Evaluated sources & licenses

Every source license/terms page was checked on `2026-11-18`, except ZINC, whose terms page was unreachable and was checked indirectly (`license_checked = "2026-11-18"`, `license_check_method = "indirect"` in `buyormake/sources.py`). Non-redistributable datasets are never committed to the repository; fetch scripts stream via `.part` temporary files with checksum verification into the gitignored `data/` directory.

| Source | Role | License | Commercial use | Redistribution | Status | Evidence |
|---|---|---|---|---|---|---|
| PubChem PUG-REST & PUG-View | Live chemical vendors & catalog IDs (no prices) | `custom/terms` (NCBI public domain + per-depositor terms) | Yes (per depositor) | No (cached locally only) | `integrated` | [NCBI Policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) |
| ZINC (ZINC15/ZINC20/ZINC22) | Purchasable / in-stock SMILES subsets | `custom/terms` | Yes ("free to use for everyone") | No (no major portions without written permission) | `fetch_only` | [ZINC Terms](https://wiki.docking.org/index.php/Terms_And_Conditions) |
| Enamine Building Blocks | Building block SDFs (manual registration) | `custom/terms` | Yes for own use after registration ("Open Structure Data" section) | No (no public republishing) | `fetch_only` | [Enamine Terms](https://enamine.net/terms-of-use) |
| eMolecules Free Download | Free SMILES + IDs (350 MB; prices require paid ePlus) | `custom/terms` | No (full commercial requires ePlus) | No | `fetch_only` | [eMolecules Readme](https://downloads.emolecules.com/free/2026-11-01/Readme.txt) |
| Mcule Database | Free in-stock & BB SMILES + IDs (no prices) | `custom/terms` | Yes (own use) | No | `fetch_only` | [Mcule Terms](https://mcule.com/terms/) |
| AiZynthFinder Code & Models | MCTS retrosynthesis + USPTO policy models | `MIT` (code) / `CC-BY-4.0` (models) | Yes | Yes | `fetch_only` | [AiZynthFinder LICENSE](https://raw.githubusercontent.com/MolecularAI/aizynthfinder/master/LICENSE) |
| AiZynthFinder ZINC Stock | `zinc_stock_17_04_20.hdf5` (663 MB, ZINC-derived) | `custom/terms` (Figshare says MIT, ZINC-derived) | No | No | `fetch_only` | [Figshare 12334577](https://figshare.com/articles/dataset/12334577) |
| USPTO Reactions (Lowe) | Reported patent reactions (1976–Sep 2016) | `CC0-1.0` | Yes | Yes | `fetch_only` | [Figshare 5104873](https://api.figshare.com/v2/articles/5104873) |
| Open Reaction Database (`ord-data`) | Reported reactions (parquet shards on Hugging Face) | `CC-BY-SA-4.0` | Yes | Yes (ShareAlike) | `fetch_only` | [ord-data LICENSE](https://raw.githubusercontent.com/open-reaction-database/ord-data/main/LICENSE) |
| RDKit Contrib SA_Score | Synthetic accessibility heuristic (1–10) | `BSD-3-Clause` | Yes | Yes | `integrated` | [RDKit license.txt](https://raw.githubusercontent.com/rdkit/rdkit/master/license.txt) |
| SCScore (Coley et al.) | Synthetic complexity score (1–5, NumPy) | `MIT` | Yes | Yes | `integrated` | [SCScore LICENSE](https://raw.githubusercontent.com/connorcoley/scscore/master/LICENSE) |
| RAscore (Reymond Group) | Learned retrosynthetic accessibility (52 MB) | `MIT` | Yes | Yes | `not_used` | [RAscore LICENSE](https://raw.githubusercontent.com/reymond-group/RAscore/master/LICENSE) |

Note on ZINC check: `wiki.docking.org` dropped TLS connections from the test host on `2026-11-18`; its terms (`fetch_only`, no major-portion redistribution without written permission) were verified from prior documentation and search results. This is machine-visible on the `Source` record via `license_check_method = "indirect"` and `notes` starting with `LICENSE CHECK INDIRECT`.

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
- **ZINC / eMolecules / Mcule / Enamine**: see [`scripts/fetch_zinc.py`](scripts/fetch_zinc.py), [`scripts/fetch_emolecules.py`](scripts/fetch_emolecules.py), [`scripts/fetch_mcule.py`](scripts/fetch_mcule.py), and [`scripts/fetch_enamine.py`](scripts/fetch_enamine.py).

## Caveats and limitations

- **Not purchasing advice.** Catalog prices change frequently and vary by region, currency, purity, pack size, institutional discounts, and account terms. A publicly listed price is not a binding quote.
- **PubChem vendor lists are capped at 25 per molecule** in the default backend (`max_offers=25`). Online cache entries default to a 30-day TTL (`--cache-max-age-days`, bypassable via `--refresh-cache`), while `--offline` uses cached entries without expiration. Transient HTTP `5xx`/`429`/timeout errors are retried and never cached as negative results.
- **Missing prices stay `null`.** Free public sources (PubChem Chemical Vendors, eMolecules free, Mcule free, ZINC, Enamine Open Structure Data) do not publish prices. Unless you supply a local offer file with real `listed_catalog` prices, `best_price` and `est_cost` remain `null`.
- **Route cost is a rough lower bound.** When leaf prices are supplied, `est_cost` sums `cheapest_USD_per_g * count * (MW_leaf / MW_target)` per 1 g target. It ignores solvents, reagents, catalysts, reaction yield, purification, analytical QC, and labor. If any leaf lacks a price, `cost_complete` is `false`. Leaves left out of the sum (unparseable SMILES or no listed USD/g price) are listed in the route's `est_cost_method` and in the `reasons` of a makeable result. For buyable molecules `est_cost` is `null` and `cost_complete` is `false`; the listed price is in `best_price`.
- **`est_cost` for makeable molecules is the cost of a found route, not the minimum-cost route.** The template search is not exhaustive: it carries at most 6 first-step disconnections into a second step, tries at most 6 second-step disconnections for each, and stops expanding a first step once one second step reaches purchasable leaves. The reported route is the best of those enumerated (fewest steps, then cost); a cheaper route may exist.
- **Prices in different currencies or pack sizes are never compared.** `best_price` is the cheapest per-gram price in `USD/g` when mass pack sizes are known; otherwise raw pack prices are compared only when all offers share one pack size, and non-USD offers only when they all use a single currency. In cross-molecule ranking and `filter_proposals(budget=...)`, `USD/g` offers are never numerically interleaved with raw per-pack `USD` or non-USD offers, and `--require-complete-cost` rejects unpriced or incomplete-cost candidates under strict budgets.
- **Template fallback is heuristic.** The 17 built-in RDKit SMARTS disconnections (15 acyclic `!@` bimolecular templates + 2 monocyclic `@` ring-openings, max 2 steps, verified for exact heavy-atom conservation) ignore stereochemistry, protecting groups, regioselectivity, chemoselectivity, and reaction conditions.
- **SA/SC scores are heuristics**, not synthesis guarantees.
- **`"unknown"` does not mean impossible.** It only means no stock offer or complete route within `max_steps` was found in the configured backends.

## License & author

- Code: MIT License ([LICENSE](LICENSE)).
- Data & model weights: each source retains its own license (see table above).
- Author: Yani Guan.
