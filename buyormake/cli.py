"""Command-line interface for buy-or-make."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

from .cache import DEFAULT_CACHE_PATH
from .rank import evaluate_candidates, filter_decisions_from_ranked
from .routes import DEFAULT_REACTION_DB
from .schema import AvailabilityRecord, export_json_schema
from .sources import SOURCES


def _read_smiles_entries(path: Path) -> tuple[list[str], list[Optional[str]]]:
    smiles_list: list[str] = []
    names_list: list[Optional[str]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split(None, 1)
        smiles_list.append(parts[0])
        names_list.append(parts[1].strip() if len(parts) > 1 and parts[1].strip() else None)
    return smiles_list, names_list


def _read_smiles_file(path: Path) -> list[str]:
    smiles_list, _ = _read_smiles_entries(path)
    return smiles_list


def _write_records_jsonl(records: list[AvailabilityRecord], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fout:
        for rec in records:
            fout.write(rec.model_dump_json() + "\n")


def _format_price(q) -> str:
    if q is None or q.value is None:
        return "unknown"
    return f"{q.value:.2f} {q.unit}"


def _add_eval_arguments(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "smiles_file",
        type=Path,
        help="Text/.smi file with one SMILES (and optional whitespace-separated molecule ID/name) per line.",
    )
    p.add_argument("--max-steps", type=int, default=3, help="Max retrosynthesis steps (default: 3).")
    p.add_argument(
        "--backends",
        default="local_stock,pubchem",
        help="Comma-separated stock backends (default: local_stock,pubchem).",
    )
    p.add_argument(
        "--route-backend",
        default="auto",
        choices=["auto", "template", "lookup", "aizynthfinder", "retro_star"],
        help="Retrosynthesis backend (default: auto).",
    )
    p.add_argument(
        "--planner",
        default="lookup_smarts",
        choices=["lookup_smarts", "retro_star"],
        help="Retrosynthesis search planner: lookup_smarts (default) or retro_star (best-first AND-OR).",
    )
    p.add_argument(
        "--reaction-db",
        type=Path,
        default=DEFAULT_REACTION_DB,
        help="Path to SQLite indexed reaction database (default: data/reactions.sqlite).",
    )
    p.add_argument("--offline", action="store_true", help="Use local cache/files only (no HTTP).")
    p.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Bypass existing cached API entries and fetch fresh vendor records.",
    )
    p.add_argument(
        "--cache-max-age-days",
        type=float,
        default=None,
        help="Maximum cache entry age in days before re-fetching in online mode (default: 30).",
    )
    p.add_argument(
        "--cache-path",
        type=Path,
        default=DEFAULT_CACHE_PATH,
        help="Path to SQLite API response cache (default: data/cache.sqlite).",
    )
    p.add_argument("--stock", action="append", default=None, help="Additional local stock file(s).")
    p.add_argument(
        "--records-out",
        type=Path,
        default=None,
        help="Write sealed AvailabilityRecord objects as JSONL to this path.",
    )
    p.add_argument(
        "--budget",
        type=float,
        default=None,
        help="Optional per-gram budget threshold for proposal filtering.",
    )
    p.add_argument(
        "--budget-currency",
        default="USD",
        help="Currency code for price comparison and --budget (default: USD).",
    )
    p.add_argument(
        "--require-complete-cost",
        action="store_true",
        help="Reject unpriced buyables and incomplete-cost routes when --budget is set.",
    )
    p.add_argument(
        "--stereo-relaxed",
        action="store_true",
        help="Allow 14-char skeletal InChIKey stereoisomer/racemate stock matches with explicit disclosure.",
    )
    p.add_argument(
        "--step-yield",
        type=float,
        default=1.0,
        help="Assumed fractional yield per reaction step in (0, 1] for route cost compounding (default: 1.0).",
    )
    p.add_argument(
        "--max-outlay-budget",
        type=float,
        default=None,
        help="Optional discrete pack MOQ cash-outlay ceiling in budget currency.",
    )
    p.add_argument("--json", action="store_true", dest="as_json", help="Output JSON.")


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="buy-or-make", description="Molecule buy-or-make triage.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_rank = sub.add_parser("rank", help="Rank SMILES from a file (with optional --budget / --records-out).")
    _add_eval_arguments(p_rank)

    p_filter = sub.add_parser("filter", help="Filter candidate proposals under availability and --budget rules.")
    _add_eval_arguments(p_filter)

    p_schema = sub.add_parser("export-schema", help="Export JSON Schema files to schema/.")
    p_schema.add_argument("--out", type=Path, default=Path("schema"), help="Output directory.")

    p_sources = sub.add_parser("sources", help="List evaluated data sources and license statuses.")
    p_sources.add_argument("--json", action="store_true", dest="as_json", help="Output JSON.")

    args = parser.parse_args(argv)

    if args.command == "export-schema":
        paths = export_json_schema(args.out)
        for p in paths:
            print(p)
        return 0

    if args.command == "sources":
        if args.as_json:
            print(json.dumps([s.model_dump(mode="json") for s in SOURCES], indent=2))
            return 0
        for s in SOURCES:
            print(
                f"{s.id:26s} status={s.status:10s} license={s.license:14s} "
                f"spdx={(s.spdx_license or '-'):24s} check={s.license_check_method:8s} "
                f"redist={str(s.redistribution_ok):5s} comm={str(s.commercial_use_ok):5s}"
            )
        return 0

    if args.command in {"rank", "filter"}:
        if not args.smiles_file.exists():
            print(f"error: file not found: {args.smiles_file}", file=sys.stderr)
            return 2
        candidates, names = _read_smiles_entries(args.smiles_file)
        backends = tuple(b.strip() for b in args.backends.split(",") if b.strip())
        ttl_sec = (
            args.cache_max_age_days * 86400.0
            if args.cache_max_age_days is not None
            else None
        )
        results, records = evaluate_candidates(
            candidates,
            candidate_names=names,
            max_steps=args.max_steps,
            backends=backends,
            route_backend=args.route_backend,
            planner=args.planner,
            offline=args.offline,
            cache_path=args.cache_path,
            cache_max_age_seconds=ttl_sec,
            refresh_cache=args.refresh_cache,
            stock_paths=args.stock,
            reaction_db_path=args.reaction_db,
            stereo_relaxed=args.stereo_relaxed,
            default_step_yield=args.step_yield,
            currency=args.budget_currency,
        )
        if args.records_out is not None:
            _write_records_jsonl(records, args.records_out)

        use_filter = (
            args.command == "filter"
            or args.budget is not None
            or args.max_outlay_budget is not None
        )
        if use_filter:
            summary = filter_decisions_from_ranked(
                results,
                budget=args.budget,
                budget_currency=args.budget_currency,
                require_complete_cost=args.require_complete_cost,
                max_outlay_budget=args.max_outlay_budget,
            )
            if args.as_json:
                print(json.dumps(summary.model_dump(mode="json"), indent=2))
            else:
                print(
                    f"{'KEEP':5s} {'STATUS':9s} {'STEPS':5s} {'PRICE':16s} "
                    f"{'EST_COST':16s} {'CANDIDATE'}"
                )
                for dec in summary.decisions:
                    r = dec.result
                    steps_s = "-" if r.n_steps is None else str(r.n_steps)
                    label = f"{r.candidate} ({r.name})" if r.name else r.candidate
                    print(
                        f"{('yes' if dec.keep else 'no'):5s} {r.status:9s} {steps_s:5s} "
                        f"{_format_price(r.best_price):16s} {_format_price(r.est_cost):16s} {label}"
                    )
            return 0

        if args.as_json:
            print(json.dumps([r.model_dump(mode="json") for r in results], indent=2))
        else:
            print(f"{'STATUS':9s} {'STEPS':5s} {'PRICE':16s} {'EST_COST':16s} {'SA':6s} {'CANDIDATE'}")
            for r in results:
                steps_s = "-" if r.n_steps is None else str(r.n_steps)
                sa_s = (
                    f"{r.sa_scores.sascore:.2f}"
                    if (r.sa_scores and r.sa_scores.sascore is not None)
                    else "-"
                )
                label = f"{r.candidate} ({r.name})" if r.name else r.candidate
                print(
                    f"{r.status:9s} {steps_s:5s} {_format_price(r.best_price):16s} "
                    f"{_format_price(r.est_cost):16s} {sa_s:6s} {label}"
                )
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
