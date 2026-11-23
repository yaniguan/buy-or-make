#!/usr/bin/env python3
"""Run buy-or-make on 20 real, novel, and invalid molecules on CPU (cold vs cached)."""

from __future__ import annotations

import argparse
import platform
import time
from pathlib import Path

import rdkit

from buyormake import RankResult, cost_aware_rank

CANDIDATES: list[tuple[str, str, str]] = [
    # (label, category, smiles)
    ("Aspirin", "common drug", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Paracetamol", "common drug", "CC(=O)Nc1ccc(O)cc1"),
    ("Ibuprofen", "common drug", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("Caffeine", "common drug", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
    ("Lidocaine", "common drug", "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    ("Benzocaine", "common drug", "CCOC(=O)c1ccc(N)cc1"),
    ("4-Bromoaniline", "building block", "Nc1ccc(Br)cc1"),
    ("Salicylic acid", "building block", "O=C(O)c1ccccc1O"),
    ("4-Methoxybenzoic acid", "building block", "COc1ccc(C(=O)O)cc1"),
    ("Pyridin-2-ylboronic acid", "building block", "OB(O)c1ccccn1"),
    ("2-Aminothiazole", "reagent/BB", "Nc1nccs1"),
    ("4-Fluorobenzenesulfonyl chloride", "reagent/BB", "O=S(=O)(Cl)c1ccc(F)cc1"),
    ("Analog-1 (fluoro-methyl benzamide)", "novel analog", "COc1cc(C)ccc1C(=O)Nc1ccc(F)c(C)c1"),
    ("Analog-2 (methoxy-methyl sulfonamide)", "novel analog", "COc1ccc(NS(=O)(=O)c2ccc(F)cc2)cc1C"),
    ("Analog-3 (biaryl amide)", "novel analog", "COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1"),
    ("Analog-4 (ester sulfonamide)", "novel analog", "CCOC(=O)c1ccc(NS(=O)(=O)c2ccc(F)cc2)cc1"),
    ("Analog-5 (cage cyclopropyl ketone)", "novel analog", "CC12CC3CC(C1)(C3)C2C(=O)C4CC4"),
    ("Analog-6 (spiro-pyrimidine lactam)", "novel analog", "O=C1NC2(CC1)CCN(c1ncccn1)CC2"),
    ("Invalid-1 (broken paren)", "invalid", "C1CC(invalid"),
    ("Invalid-2 (unclosed ring)", "invalid", "c1ccccc"),
]


def _cpu_model() -> str:
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        for line in cpuinfo.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    return platform.processor() or platform.machine()


def _fmt_price(q) -> str:
    if q is None or q.value is None:
        return "null"
    return f"{q.value:.2f} {q.unit}"


def _coverage_stats(results: list[RankResult]) -> dict[str, int]:
    total = len(results)
    buyable = sum(1 for r in results if r.status == "buyable")
    makeable = sum(1 for r in results if r.status == "makeable")
    unknown = sum(1 for r in results if r.status == "unknown")
    invalid = sum(1 for r in results if r.status == "invalid")
    valid = total - invalid
    price_known = sum(
        1
        for r in results
        if (r.best_price is not None and r.best_price.value is not None)
        or (r.est_cost is not None and r.est_cost.value is not None and r.cost_complete)
    )
    route_found = makeable
    return {
        "total": total,
        "valid": valid,
        "buyable": buyable,
        "makeable": makeable,
        "unknown": unknown,
        "invalid": invalid,
        "price_known": price_known,
        "route_found": route_found,
    }


def render_markdown(
    results: list[RankResult],
    cold_sec: float,
    cached_sec: float,
    cpu: str,
    offline_only: bool,
) -> str:
    label_by_smiles = {smi: (lbl, cat) for lbl, cat, smi in CANDIDATES}
    stats = _coverage_stats(results)
    n = len(CANDIDATES)
    mode_str = "offline cache only" if offline_only else "live PubChem PUG-REST/PUG-View + local SA/SCScore + template fallback"

    lines = [
        "# buy-or-make Demo Run Output",
        "",
        f"- **Mode**: {mode_str}",
        f"- **CPU**: {cpu}",
        f"- **Python / RDKit**: Python {platform.python_version()} / RDKit {rdkit.__version__}",
        f"- **Candidates**: {n} ({stats['valid']} valid, {stats['invalid']} invalid)",
        f"- **Cold wall time**: {cold_sec:.3f} s ({cold_sec / n:.3f} s / molecule)",
        f"- **Cached wall time**: {cached_sec:.3f} s ({cached_sec / n:.3f} s / molecule)",
        "",
        "## Coverage Summary",
        "",
        "| Metric | Count | Fraction of total | Fraction of valid |",
        "|---|---:|---:|---:|",
    ]
    for key, label in [
        ("buyable", "Buyable (listed by >=1 vendor)"),
        ("makeable", "Makeable (complete route <= 2 steps)"),
        ("unknown", "Unknown (not in stock; no complete route <= 2 steps)"),
        ("invalid", "Invalid SMILES"),
        ("price_known", "Complete price known"),
        ("route_found", "Route found (for non-buyable)"),
    ]:
        c = stats[key]
        frac_tot = f"{(100.0 * c / stats['total']):.1f}%"
        frac_val = f"{(100.0 * c / stats['valid']):.1f}%" if key != "invalid" and stats["valid"] else "-"
        lines.append(f"| {label} | {c} | {frac_tot} | {frac_val} |")

    lines.extend(
        [
            "",
            "## Ranked Candidates",
            "",
            "| # | Label | Category | Status | Steps | Best Price | Est. Cost | SA | SC | Evidence URLs | Notes |",
            "|---:|---|---|---|---:|---|---|---:|---:|---:|---|",
        ]
    )
    for idx, r in enumerate(results, start=1):
        lbl, cat = label_by_smiles.get(r.candidate, (r.candidate, ""))
        steps_s = "-" if r.n_steps is None else str(r.n_steps)
        sa_s = f"{r.sa_scores.sascore:.2f}" if (r.sa_scores and r.sa_scores.sascore is not None) else "-"
        sc_s = f"{r.sa_scores.scscore:.2f}" if (r.sa_scores and r.sa_scores.scscore is not None) else "-"
        note = r.reasons[0] if r.reasons else ""
        lines.append(
            f"| {idx} | {lbl} | {cat} | `{r.status}` | {steps_s} | "
            f"{_fmt_price(r.best_price)} | {_fmt_price(r.est_cost)} | {sa_s} | {sc_s} | "
            f"{len(r.evidence)} | {note} |"
        )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run 20-molecule buy-or-make CPU demo.")
    parser.add_argument("--cache", type=Path, default=Path("data/demo_cache.sqlite"))
    parser.add_argument("--scscore-weights", type=Path, default=None)
    parser.add_argument("--save-md", type=Path, default=Path("examples/demo_output.md"))
    parser.add_argument("--offline", action="store_true", help="Skip live PubChem cold run.")
    args = parser.parse_args()

    smiles_list = [smi for _, _, smi in CANDIDATES]
    cpu = _cpu_model()
    extra_kw = {"scscore_weights": args.scscore_weights} if args.scscore_weights else {}

    if not args.offline and args.cache.exists():
        args.cache.unlink()

    t0 = time.perf_counter()
    cost_aware_rank(
        smiles_list,
        max_steps=2,
        backends=("local_stock", "pubchem"),
        route_backend="auto",
        offline=args.offline,
        cache_path=args.cache,
        stock_paths=[],
        **extra_kw,
    )
    cold_sec = time.perf_counter() - t0

    t1 = time.perf_counter()
    results_cached = cost_aware_rank(
        smiles_list,
        max_steps=2,
        backends=("local_stock", "pubchem"),
        route_backend="auto",
        offline=True,
        cache_path=args.cache,
        stock_paths=[],
        **extra_kw,
    )
    cached_sec = time.perf_counter() - t1

    md = render_markdown(results_cached, cold_sec, cached_sec, cpu, offline_only=args.offline)
    print(md)
    if args.save_md:
        args.save_md.parent.mkdir(parents=True, exist_ok=True)
        args.save_md.write_text(md, encoding="utf-8")
        print(f"Saved report to {args.save_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
