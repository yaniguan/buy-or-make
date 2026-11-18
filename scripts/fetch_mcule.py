#!/usr/bin/env python3
"""Fetch Mcule free in-stock building blocks SMILES into data/mcule/ (fetch_only; no prices)."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
import httpx

PAGE_URL = "https://mcule.com/database/"
OUT_DIR = Path("data/mcule")


def find_in_stock_bb_url() -> str | None:
    r = httpx.get(PAGE_URL, follow_redirects=True, timeout=20.0)
    r.raise_for_status()
    m = re.search(
        r'https://[^\s"\']*mcule_purchasable_in_stock_building_blocks_[^\s"\']*\.smi\.gz',
        r.text,
    )
    if m:
        return m.group(0)
    m2 = re.search(r'href="([^"]*mcule_purchasable_in_stock_building_blocks_[^"]*\.smi\.gz)"', r.text)
    if m2:
        href = m2.group(1)
        return href if href.startswith("http") else f"https://mcule.com{href}"
    return None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch Mcule in-stock building blocks SMILES (~151 MB; structures + Mcule IDs only)."
    )
    parser.add_argument("--dry-run", action="store_true", help="Print resolved URL and exit.")
    args = parser.parse_args()

    print(
        "Mcule notice: free database downloads include 2D structures and Mcule IDs only (no prices). "
        "Terms allow own business/research use only; no third-party redistribution (mcule.com/terms/)."
    )
    url = find_in_stock_bb_url()
    if not url:
        print("Could not auto-detect download link on https://mcule.com/database/ (download manually into data/mcule/).")
        return 1
    if args.dry_run:
        print(f"Resolved URL: {url}")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dest = OUT_DIR / Path(url).name
    print(f"Downloading {url} -> {dest}")
    with httpx.stream("GET", url, follow_redirects=True, timeout=300.0) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_bytes():
                f.write(chunk)
    print(f"Saved {dest} ({dest.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
