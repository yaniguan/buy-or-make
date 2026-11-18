# buy-or-make Demo Run Output

- **Mode**: live PubChem PUG-REST/PUG-View + local SA/SCScore + template fallback
- **CPU**: Intel(R) Xeon(R) Processor
- **Python / RDKit**: Python 3.12.16 / RDKit 2026.09.2
- **Candidates**: 20 (18 valid, 2 invalid)
- **Cold wall time**: 24.970 s (1.249 s / molecule)
- **Cached wall time**: 0.099 s (0.005 s / molecule)

## Coverage Summary

| Metric | Count | Fraction of total | Fraction of valid |
|---|---:|---:|---:|
| Buyable (listed by >=1 vendor) | 14 | 70.0% | 77.8% |
| Makeable (complete route <= 2 steps) | 3 | 15.0% | 16.7% |
| Unknown (not in stock; no complete route <= 2 steps) | 1 | 5.0% | 5.6% |
| Invalid SMILES | 2 | 10.0% | - |
| Complete price known | 0 | 0.0% | 0.0% |
| Route found (for non-buyable) | 3 | 15.0% | 16.7% |

## Ranked Candidates

| # | Label | Category | Status | Steps | Best Price | Est. Cost | SA | SC | Evidence URLs | Notes |
|---:|---|---|---|---:|---|---|---:|---:|---:|---|
| 1 | Aspirin | common drug | `buyable` | 0 | null | null | 1.58 | 1.59 | 25 | listed by 25 catalog offer(s): CATO Research Chemicals Inc. (+22 more offers) |
| 2 | Paracetamol | common drug | `buyable` | 0 | null | null | 1.41 | 1.54 | 25 | listed by 25 catalog offer(s): LGC Standards (+22 more offers) |
| 3 | Ibuprofen | common drug | `buyable` | 0 | null | null | 2.19 | 2.60 | 25 | listed by 25 catalog offer(s): AppChem, A2B Chem (+22 more offers) |
| 4 | Caffeine | common drug | `buyable` | 0 | null | null | 2.30 | 2.79 | 25 | listed by 25 catalog offer(s): Thermo Fisher Scientific (+22 more offers) |
| 5 | Lidocaine | common drug | `buyable` | 0 | null | null | 1.81 | 2.26 | 24 | listed by 25 catalog offer(s): 001Chemical, Accela ChemBio Inc., MolCore (+22 more offers) |
| 6 | Benzocaine | common drug | `buyable` | 0 | null | null | 1.45 | 1.50 | 25 | listed by 25 catalog offer(s): Thermo Fisher Scientific (+22 more offers) |
| 7 | 4-Bromoaniline | building block | `buyable` | 0 | null | null | 1.54 | 1.19 | 25 | listed by 25 catalog offer(s): Enamine, abcr GmbH (+22 more offers) |
| 8 | Salicylic acid | building block | `buyable` | 0 | null | null | 1.43 | 1.46 | 25 | listed by 25 catalog offer(s): LGC Standards (+22 more offers) |
| 9 | 4-Methoxybenzoic acid | building block | `buyable` | 0 | null | null | 1.17 | 1.56 | 25 | listed by 25 catalog offer(s): ChemFaces, AA BLOCKS, SynQuest Laboratories (+22 more offers) |
| 10 | Pyridin-2-ylboronic acid | building block | `buyable` | 0 | null | null | 2.96 | 1.46 | 25 | listed by 25 catalog offer(s): Finetech Industry Limited, BOC Sciences (+22 more offers) |
| 11 | 2-Aminothiazole | reagent/BB | `buyable` | 0 | null | null | 3.17 | 1.49 | 25 | listed by 25 catalog offer(s): Accela ChemBio Inc., Enamine (+22 more offers) |
| 12 | 4-Fluorobenzenesulfonyl chloride | reagent/BB | `buyable` | 0 | null | null | 1.68 | 1.19 | 25 | listed by 25 catalog offer(s): Thermo Fisher Scientific (+22 more offers) |
| 13 | Analog-2 (methoxy-methyl sulfonamide) | novel analog | `buyable` | 0 | null | null | 1.64 | 2.31 | 1 | listed by 1 catalog offer(s): AKos Consulting & Solutions |
| 14 | Analog-4 (ester sulfonamide) | novel analog | `buyable` | 0 | null | null | 1.62 | 2.02 | 25 | listed by 25 catalog offer(s): A2B Chem, CHIRALEN, AppChem (+22 more offers) |
| 15 | Analog-1 (fluoro-methyl benzamide) | novel analog | `makeable` | 1 | null | null | 1.66 | 2.36 | 50 | makeable in 1 step(s) via template_fallback (heuristic RDKit SMARTS disconnections, max 2 steps; ignores stereochemistry and selectivity) from 2 purchasable leaf/leaves |
| 16 | Analog-3 (biaryl amide) | novel analog | `makeable` | 1 | null | null | 1.54 | 3.21 | 49 | makeable in 1 step(s) via template_fallback (heuristic RDKit SMARTS disconnections, max 2 steps; ignores stereochemistry and selectivity) from 2 purchasable leaf/leaves |
| 17 | Analog-6 (spiro-pyrimidine lactam) | novel analog | `makeable` | 1 | null | null | 3.17 | 3.47 | 27 | makeable in 1 step(s) via template_fallback (heuristic RDKit SMARTS disconnections, max 2 steps; ignores stereochemistry and selectivity) from 2 purchasable leaf/leaves |
| 18 | Analog-5 (cage cyclopropyl ketone) | novel analog | `unknown` | - | null | null | 5.49 | 3.03 | 0 | not found in stock backends and no complete route to purchasable leaves found within max_steps=2 (unknown does not mean unmakeable) |
| 19 | Invalid-1 (broken paren) | invalid | `invalid` | - | null | null | - | - | 0 | invalid SMILES; RDKit could not parse structure |
| 20 | Invalid-2 (unclosed ring) | invalid | `invalid` | - | null | null | - | - | 0 | invalid SMILES; RDKit could not parse structure |
