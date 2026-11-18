"""Command-line interface for buy-or-make."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

from .rank import cost_aware_rank
from .schema import export_json_schema
from .sources import SOURCES


def _read_smiles_file(path: Path) -> list[str]:
    lines: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        lines.append(s.split()[0])
    return lines


def _format_price(q) -> str:
    if q is None or q.value is None:
        return "unknown"
    return f"{q.value:.2f} {q.unit}"


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="buy-or-make", description="Molecule buy-or-make triage.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_rank = sub.add_parser("rank", help="Rank SMILES from a file.")
    p_rank.add_argument("smiles_file", type=Path, help="Text file with one SMILES per line.")
    p_rank.add_argument("--max-steps", type=int, default=3, help="Max retrosynthesis steps (default: 3).")
    p_rank.add_argument(
        "--backends",
        default="local_stock,pubchem",
        help="Comma-separated stock backends (default: local_stock,pubchem).",
    )
    p_rank.add_argument(
        "--route-backend",
        default="auto",
        choices=["auto", "template", "lookup", "aizynthfinder"],
        help="Retrosynthesis backend (default: auto).",
    )
    p_rank.add_argument("--offline", action="store_true", help="Use local cache/files only (no HTTP).")
    p_rank.add_argument("--stock", action="append", default=None, help="Additional local stock file(s).")
    p_rank.add_argument("--json", action="store_true", dest="as_json", help="Output JSON.")

    p_schema = sub.add_parser("export-schema", help="Export JSON Schema files to schema/.")
    p_schema.add_argument("--out", type=Path, default=Path("schema"), help="Output directory.")

    sub.add_parser("sources", help="List evaluated data sources and license statuses.")

    args = parser.parse_args(argv)

    if args.command == "export-schema":
        paths = export_json_schema(args.out)
        for p in paths:
            print(p)
        return 0

    if args.command == "sources":
        for s in SOURCES:
            print(
                f"{s.id:26s} status={s.status:10s} license={s.license:14s} "
                f"redist={str(s.redistribution_ok):5s} comm={str(s.commercial_use_ok):5s}"
            )
        return 0

    if args.command == "rank":
        if not args.smiles_file.exists():
            print(f"error: file not found: {args.smiles_file}", file=sys.stderr)
            return 2
        candidates = _read_smiles_file(args.smiles_file)
        backends = tuple(b.strip() for b in args.backends.split(",") if b.strip())
        results = cost_aware_rank(
            candidates,
            max_steps=args.max_steps,
            backends=backends,
            route_backend=args.route_backend,
            offline=args.offline,
            stock_paths=args.stock,
        )
        if args.as_json:
            print(json.dumps([r.model_dump(mode="json") for r in results], indent=2))
        else:
            print(f"{'STATUS':9s} {'STEPS':5s} {'PRICE':16s} {'EST_COST':16s} {'SA':6s} {'CANDIDATE'}")
            for r in results:
                steps_s = "-" if r.n_steps is None else str(r.n_steps)
                sa_s = f"{r.sa_scores.sascore:.2f}" if (r.sa_scores and r.sa_scores.sascore is not None) else "-"
                print(
                    f"{r.status:9s} {steps_s:5s} {_format_price(r.best_price):16s} "
                    f"{_format_price(r.est_cost):16s} {sa_s:6s} {r.candidate}"
                )
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
