#!/usr/bin/env python3
"""Fetch eMolecules free SMILES dump (structures + IDs only, no prices) into data/emolecules/."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
import httpx

INDEX_URL = "https://downloads.emolecules.com/free/"
OUT_DIR = Path("data/emolecules")


def latest_release_url() -> str:
    r = httpx.get(INDEX_URL, timeout=20.0)
    r.raise_for_status()
    dirs = sorted(re.findall(r'href="(\d{4}-\d{2}-\d{2}/)"', r.text))
    if not dirs:
        raise RuntimeError(f"Could not find dated release directories at {INDEX_URL}")
    return f"{INDEX_URL}{dirs[-1]}version.smi.gz"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch eMolecules free version.smi.gz (~350 MB; structures only, no prices)."
    )
    parser.add_argument("--dry-run", action="store_true", help="Print resolved release URL and exit.")
    args = parser.parse_args()

    url = latest_release_url()
    print("eMolecules notice: free download contains SMILES and IDs only (no prices, no redistribution).")
    if args.dry_run:
        print(f"Resolved URL: {url}")
        return 0

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dest = OUT_DIR / "version.smi.gz"
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
