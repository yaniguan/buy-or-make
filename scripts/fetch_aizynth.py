#!/usr/bin/env python3
"""Fetch AiZynthFinder public USPTO policy models (CC-BY-4.0) and optional ZINC stock (fetch_only)."""

from __future__ import annotations

import argparse
from pathlib import Path
import httpx

# Figshare 12334577 files (also reachable via aizynthfinder's `download_public_data`)
FILES = {
    "uspto_model.onnx": "https://zenodo.org/records/7797465/files/uspto_model.onnx?download=1",
    "uspto_templates.csv.gz": "https://zenodo.org/records/7341155/files/uspto_templates.csv.gz?download=1",
    "uspto_filter_model.onnx": "https://zenodo.org/records/7797465/files/uspto_filter_model.onnx?download=1",
    "zinc_stock_17_04_20.hdf5": "https://ndownloader.figshare.com/files/23086469",
}
OUT_DIR = Path("data/aizynth")


def main() -> int:
    parser = argparse.ArgumentParser(description="Fetch AiZynthFinder public data (~790 MB total).")
    parser.add_argument("--with-stock", action="store_true", help="Also download 663 MB ZINC HDF5 stock.")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = ["uspto_model.onnx", "uspto_templates.csv.gz", "uspto_filter_model.onnx"]
    if args.with_stock:
        targets.append("zinc_stock_17_04_20.hdf5")

    print("Note: If Zenodo returns 403 from automated clients, run `download_public_data data/aizynth` "
          "after `pip install -e .[aizynthfinder]`.")

    for name in targets:
        dest = OUT_DIR / name
        if dest.exists() and dest.stat().st_size > 100_000:
            print(f"Already present: {dest}")
            continue
        url = FILES[name]
        print(f"Downloading {url} -> {dest}")
        try:
            with httpx.stream("GET", url, follow_redirects=True, timeout=120.0) as r:
                r.raise_for_status()
                with dest.open("wb") as f:
                    for chunk in r.iter_bytes():
                        f.write(chunk)
        except httpx.HTTPError as exc:
            print(f"Failed to download {name}: {exc}")
            return 1

    cfg = OUT_DIR / "config.yml"
    cfg.write_text(
        "expansion:\n"
        "  uspto:\n"
        f"    - {OUT_DIR.resolve() / 'uspto_model.onnx'}\n"
        f"    - {OUT_DIR.resolve() / 'uspto_templates.csv.gz'}\n"
        "filter:\n"
        f"  uspto: {OUT_DIR.resolve() / 'uspto_filter_model.onnx'}\n"
        "stock:\n"
        f"  zinc: {OUT_DIR.resolve() / 'zinc_stock_17_04_20.hdf5'}\n",
        encoding="utf-8",
    )
    print(f"Wrote {cfg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
