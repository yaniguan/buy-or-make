#!/usr/bin/env python3
"""Fetch ZINC20/ZINC22 in-stock SMILES subsets into data/zinc/ (fetch_only; no redistribution)."""

from __future__ import annotations

import argparse
from pathlib import Path
import httpx

OUT_DIR = Path("data/zinc")

# ZINC20 in-stock tranche catalog index (academic/research use; see wiki.docking.org Terms_And_Conditions)
DEFAULT_TRANCHE_URLS = [
    "https://files.docking.org/2D/AA/AAAA.smi",
]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch ZINC SMILES tranches into data/zinc/ (not redistributed in repo)."
    )
    parser.add_argument(
        "--url",
        action="append",
        default=None,
        help="ZINC .smi/.smi.gz tranche URL(s) to fetch.",
    )
    args = parser.parse_args()

    print(
        "ZINC license notice: free for academic/research use; redistribution of major portions "
        "without written permission is prohibited (wiki.docking.org/index.php/Terms_And_Conditions)."
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    urls = args.url or DEFAULT_TRANCHE_URLS
    for url in urls:
        dest = OUT_DIR / Path(url).name
        if dest.exists() and dest.stat().st_size > 0:
            print(f"Already present: {dest}")
            continue
        print(f"Downloading {url} -> {dest}")
        try:
            with httpx.stream("GET", url, follow_redirects=True, timeout=60.0) as r:
                r.raise_for_status()
                with dest.open("wb") as f:
                    for chunk in r.iter_bytes():
                        f.write(chunk)
        except httpx.HTTPError as exc:
            print(f"Download failed for {url}: {exc}")
            return 1
    print(f"Done. Pass --stock {OUT_DIR / Path(urls[0]).name} to buy-or-make rank.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
