#!/usr/bin/env python3
"""Fetch ZINC20/ZINC22 in-stock SMILES subsets into data/zinc/ (fetch_only; no redistribution)."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import httpx

from buyormake.canonical import atomic_download

OUT_DIR = Path("data/zinc")

# ZINC20 in-stock tranche catalog index (academic/research use; see wiki.docking.org Terms_And_Conditions)
DEFAULT_TRANCHE_URLS = [
    "https://files.docking.org/2D/AA/AAAA.smi",
]


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Fetch ZINC SMILES tranches into data/zinc/ (not redistributed in repo)."
    )
    parser.add_argument(
        "--url",
        action="append",
        default=None,
        help="ZINC .smi/.smi.gz tranche URL(s) to fetch.",
    )
    parser.add_argument("--expected-sha256", default=None, help="Optional pinned SHA-256 hex digest.")
    parser.add_argument("--force", action="store_true", help="Re-download even if present.")
    args = parser.parse_args(argv)

    print(
        "ZINC license notice: free for academic/research use; redistribution of major portions "
        "without written permission is prohibited (wiki.docking.org/index.php/Terms_And_Conditions)."
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    urls = args.url or DEFAULT_TRANCHE_URLS
    for url in urls:
        dest = OUT_DIR / Path(url).name
        print(f"Downloading {url} -> {dest}")
        try:
            sha = atomic_download(
                url,
                dest,
                expected_sha256=args.expected_sha256 if len(urls) == 1 else None,
                timeout=60.0,
                force=args.force,
            )
            print(f"Saved {dest} ({dest.stat().st_size} bytes, sha256={sha})")
        except (httpx.HTTPError, ValueError) as exc:
            print(f"Download failed for {url}: {exc}")
            return 1
    print(f"Done. Pass --stock {OUT_DIR / Path(urls[0]).name} to buy-or-make rank.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
