# buy-or-make

Check whether candidate molecules are directly purchasable from chemical vendors or makeable in N steps from purchasable building blocks, with conservative leaf-only cost lower bounds computed strictly from real listed catalog prices. Missing prices stay `null` and unverified molecules stay `"unknown"`.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"

# Optional: fetch MIT-licensed SCScore NumPy weights (6.18 MB) into data/scscore/
python scripts/fetch_scscore.py

# Rank candidates from a SMILES file (1 SMILES per line)
cat > /tmp/candidates.smi << 'EOF'
CC(=O)Oc1ccccc1C(=O)O
COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1
invalid_smiles
EOF

buy-or-make rank /tmp/candidates.smi --max-steps 2
buy-or-make rank /tmp/candidates.smi --max-steps 2 --json
```

## Python API

```python
from buyormake import cost_aware_rank, filter_proposals

candidates = [
    "CC(=O)Oc1ccccc1C(=O)O",                    # buyable
    "COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1",    # novel analog
    "C1CC(invalid",                              # invalid
]

ranked = cost_aware_rank(
    candidates,
    max_steps=2,
    backends=("local_stock", "pubchem"),
    route_backend="auto",
    offline=False,
)
for r in ranked:
    print(r.status, r.n_steps, r.best_price, r.est_cost, r.reasons[0])

# Co-scientist proposal filter (keeps buyable/makeable, flags over-budget)
summary = filter_proposals(candidates, max_steps=2, budget=25.0)
print("kept:", len(summary.kept), "over_budget:", len(summary.over_budget))
```

Ranking order: `buyable` (cheaper known price first, `null` price after) < `makeable` (fewer steps, lower cost) < `unknown` < `invalid`.

## Schema summary

All models in [`buyormake/schema.py`](buyormake/schema.py) use Pydantic v2 (`extra="forbid"`), align shared field names (`Source`, `Quantity`, `Uncertainty`, `Provenance`, `Molecule`, `origin`, `schema_version`, `record_hash`) with the `assayfold` hub schema, and export JSON Schemas to `schema/`:

- `Source`: evaluated dataset/API metadata, SPDX or `custom/terms`, `license_evidence`, `license_quote`, `license_checked`, `commercial_use_ok`, `redistribution_ok`, `status`.
- `Offer`: `inchikey`, `vendor`, `catalog_id`, `pack_size` (`Quantity`), `purity`, `price` (`float | null`), `currency` (`str | null`), `price_basis` (`"listed_catalog" | "unknown"`), `lead_time_days`, `in_stock`, `source_id`, `source_url`, `fetched_utc`, `license`.
- `Route` / `RouteStep`: `target_inchikey`, `n_steps`, `steps`, `leaves`, `leaves_purchasable`, `method`, `est_cost` (`Quantity | null`), `est_cost_method`, `est_cost_complete`.
- `AvailabilityRecord`: `molecule`, `offers`, `routes`, `sa_scores` (`sascore`, `scscore`, `rascore`), `origin`, `provenance`, `licenses`, `schema_version`, deterministic `record_hash` (SHA-256 of canonical JSON).
- `RankResult`: `candidate`, `canonical_smiles`, `inchikey`, `status` (`buyable | makeable | unknown | invalid`), `n_steps`, `best_price`, `est_cost`, `cost_complete`, `sa_scores`, `reasons`, `evidence`.

Regenerate JSON Schemas with:

```bash
buy-or-make export-schema --out schema
```

## Evaluated sources & licenses

Every source license/terms page was checked on `2026-11-18`, except ZINC, whose terms page was unreachable and was checked indirectly (see note below; recorded as `license_checked` = `2026-11-18 (indirect: ...)` in `buyormake/sources.py`). Non-redistributable datasets are never committed to the repository; fetch scripts write into the gitignored `data/` directory.

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

Note on ZINC check: `wiki.docking.org` dropped TLS connections from the test host on `2026-11-18`; its terms (`fetch_only`, no major-portion redistribution without written permission) were verified from prior documentation and search results. This is marked in the `Source` record itself: ZINC's `license_checked` carries an `(indirect: ...)` suffix and its `notes` start with `LICENSE CHECK INDIRECT`.

## Demo coverage and CPU runtime

Measured from an actual run of [`examples/demo.py`](examples/demo.py) on `2026-11-18` (20 candidates: 6 common drugs, 6 reagents/building blocks, 6 hand-drawn analogs, 2 invalid SMILES; 2 of the 6 analogs turned out to be already listed by vendors) using live PubChem PUG-REST/PUG-View + RDKit SAscore + NumPy SCScore + built-in 2-step template fallback:

- **Hardware / Environment**: Intel(R) Xeon(R) Processor (CPU only), Python 3.12.16, RDKit 2026.09.2
- **Cold run wall time**: `24.970 s` total (`1.249 s / molecule`, rate-limited at <= 5 req/s to PubChem)
- **Cached run wall time**: `0.099 s` total (`0.005 s / molecule`, offline SQLite cache)

| Metric | Count | Fraction of total (n=20) | Fraction of valid (n=18) |
|---|---:|---:|---:|
| Buyable (listed by >=1 vendor in PubChem) | 14 | 70.0% | 77.8% |
| Makeable (complete route <= 2 steps) | 3 | 15.0% | 16.7% |
| Unknown (not in stock; no complete route <= 2 steps) | 1 | 5.0% | 5.6% |
| Invalid SMILES | 2 | 10.0% | - |
| Complete price known | 0 | 0.0% | 0.0% |
| Route found (for non-buyable) | 3 | 15.0% | 16.7% |

The demo configuration has no price source: PubChem Chemical Vendors publishes vendors and catalog IDs but not prices, and no local priced offer file was supplied. The cost columns (`best_price`, `est_cost`) are therefore structurally `null`, and the "Complete price known = 0" row is a consequence of that configuration, not a measured result of the cost logic. The cost estimator (`cheapest_offer_for_target`, `estimate_route_cost`) is demonstrated only by the unit tests in `tests/test_cost.py` and `tests/test_rank.py`, which use synthetic local offers.

Full candidate-by-candidate table is saved in [`examples/demo_output.md`](examples/demo_output.md).

## Optional backends & local data scripts

- **Local priced/unpriced stock** (`.tsv`, `.csv`, `.jsonl`, `.smi`, `.smi.gz`, `.sqlite`, `.hdf5`): pass `--stock path/to/file` or place `data/stock.sqlite` / `data/stock.tsv`.
- **USPTO / ORD exact reaction lookup**:
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
- **PubChem vendor lists are capped at 25 per molecule** in the default backend (`max_offers=25`). A PubChem listing means a vendor deposited the record; it does not prove the item is in stock today.
- **Missing prices stay `null`.** Free public sources (PubChem Chemical Vendors, eMolecules free, Mcule free, ZINC, Enamine Open Structure Data) do not publish prices. Unless you supply a local offer file with real `listed_catalog` prices, `best_price` and `est_cost` remain `null`.
- **Route cost is a rough lower bound.** When leaf prices are supplied, `est_cost` sums `cheapest_USD_per_g * count * (MW_leaf / MW_target)` per 1 g target. It ignores solvents, reagents, catalysts, reaction yield, purification, analytical QC, and labor. If any leaf lacks a price, `cost_complete` is `false`. Leaves left out of the sum (unparseable SMILES or no listed USD/g price) are listed in the route's `est_cost_method` and in the `reasons` of a makeable result. For buyable molecules `est_cost` is `null` and `cost_complete` is `false`; the listed price is in `best_price`.
- **`est_cost` for makeable molecules is the cost of a found route, not the minimum-cost route.** The template search is not exhaustive: it carries at most 6 first-step disconnections into a second step, tries at most 6 second-step disconnections for each, and stops expanding a first step once one second step reaches purchasable leaves. The reported route is the best of those enumerated (fewest steps, then cost); a cheaper route may exist.
- **Prices in different currencies or pack sizes are never compared.** `best_price` is the cheapest per-gram price in USD when mass pack sizes are known; otherwise raw pack prices are compared only when all offers share one pack size, and non-USD offers only when they all use a single currency. Anything else gives `null`.
- **Template fallback is heuristic.** The 14 built-in RDKit SMARTS disconnections (max 2 steps) ignore stereochemistry, protecting groups, regioselectivity, chemoselectivity, and reaction conditions.
- **SA/SC scores are heuristics**, not synthesis guarantees.
- **`"unknown"` does not mean impossible.** It only means no stock offer or complete route within `max_steps` was found in the configured backends.

## License & author

- Code: MIT License ([LICENSE](LICENSE)).
- Data & model weights: each source retains its own license (see table above).
- Author: Yani Guan.
