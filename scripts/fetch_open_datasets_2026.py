#!/usr/bin/env python3
"""On-demand fetcher & SQLite indexer for 2024-2026 open chemistry corpora.

All downloaded dumps stay inside gitignored `data/open_corpora/<source_id>/` and are verified
against pinned SHA-256 / MD5 / byte-length metadata before atomic replacement.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

try:
    from buyormake.open_corpora import (
        CORPUS_ARTIFACTS,
        DEFAULT_OPEN_CORPORA_DIR,
        DEFAULT_REACTION_DB,
        DEFAULT_REFERENCE_COMPOUNDS_DB,
        DEFAULT_STOCK_DB,
        download_corpus_artifact,
        ingest_downloaded_artifact,
    )
    from buyormake.sources import SOURCE_BY_ID
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from buyormake.open_corpora import (
        CORPUS_ARTIFACTS,
        DEFAULT_OPEN_CORPORA_DIR,
        DEFAULT_REACTION_DB,
        DEFAULT_REFERENCE_COMPOUNDS_DB,
        DEFAULT_STOCK_DB,
        download_corpus_artifact,
        ingest_downloaded_artifact,
    )
    from buyormake.sources import SOURCE_BY_ID


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Fetch and index audited 2024-2026 open stock, compound, and retrosynthesis corpora."
    )
    parser.add_argument(
        "--list",
        action="store_true",
        dest="list_only",
        help="List all registered open-corpus artifacts with SPDX licenses and sizes.",
    )
    parser.add_argument(
        "--artifact",
        choices=sorted(CORPUS_ARTIFACTS.keys()),
        default=None,
        help="Registered artifact key to download and optionally index.",
    )
    parser.add_argument(
        "--local-file",
        type=Path,
        default=None,
        help="Ingest an already-downloaded local file for `--artifact` without network I/O.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OPEN_CORPORA_DIR,
        help="Root directory for downloaded dumps (default: data/open_corpora).",
    )
    parser.add_argument(
        "--stock-db",
        type=Path,
        default=DEFAULT_STOCK_DB,
        help="Destination SQLite for stopping-stock building blocks (default: data/stock.sqlite).",
    )
    parser.add_argument(
        "--reference-db",
        type=Path,
        default=DEFAULT_REFERENCE_COMPOUNDS_DB,
        help="Destination SQLite for reference compound corpora (default: data/reference_compounds.sqlite).",
    )
    parser.add_argument(
        "--reaction-db",
        type=Path,
        default=DEFAULT_REACTION_DB,
        help="Destination SQLite for reaction/route corpora (default: data/reactions.sqlite).",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=50_000,
        help="Max rows to import/index (0 = all; default: 50000).",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Download and verify checksum only; skip SQLite indexing.",
    )
    parser.add_argument(
        "--expected-sha256",
        default=None,
        help="Optional explicit SHA-256 hex digest override.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download artifact even if local file already passes checksum verification.",
    )
    args = parser.parse_args(argv)

    if args.list_only or args.artifact is None:
        for key in sorted(CORPUS_ARTIFACTS):
            art = CORPUS_ARTIFACTS[key]
            src = SOURCE_BY_ID[art.source_id]
            size_mb = f"{art.expected_size / 1_000_000:.2f} MB" if art.expected_size else "variable"
            print(
                f"{key:32s} kind={art.kind:19s} spdx={(src.spdx_license or src.license):12s} "
                f"size={size_mb:12s} {art.description}"
            )
        return 0

    if args.local_file is not None:
        if not args.local_file.exists():
            print(f"File not found: {args.local_file}")
            return 1
        target_path = args.local_file
    else:
        target_path = download_corpus_artifact(
            args.artifact,
            out_root=args.out_dir,
            expected_sha256=args.expected_sha256,
            force=args.force,
        )
        print(f"Verified download: {target_path}")

    if args.download_only:
        return 0

    count = ingest_downloaded_artifact(
        args.artifact,
        target_path,
        stock_db=args.stock_db,
        reference_db=args.reference_db,
        reaction_db=args.reaction_db,
        max_rows=args.max_rows,
    )
    print(f"Ingested {count} records for {args.artifact} from {target_path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
