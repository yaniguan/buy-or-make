#!/usr/bin/env python3
"""Fetch AiZynthFinder public USPTO policy models (CC-BY-4.0) and optional ZINC stock (fetch_only)."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import httpx

from buyormake.canonical import atomic_download

# Zenodo 7797465 / 7341155 & Figshare 12334577 official metadata (URL, byte size, MD5)
FILES: dict[str, dict[str, object]] = {
    "uspto_model.onnx": {
        "url": "https://zenodo.org/records/7797465/files/uspto_model.onnx?download=1",
        "size": 91_518_243,
        "md5": "b466a436b6a8eceaba8c14ef5e2a88c8",
    },
    "uspto_templates.csv.gz": {
        "url": "https://zenodo.org/records/7341155/files/uspto_unique_templates.csv.gz?download=1",
        "size": 3_313_598,
        "md5": "e8fb29d4dcd54073658000e212b7845a",
    },
    "uspto_filter_model.onnx": {
        "url": "https://zenodo.org/records/7797465/files/uspto_filter_model.onnx?download=1",
        "size": 16_788_859,
        "md5": "ab0df54b62a0aeceecd68051f8dc3f2e",
    },
    "zinc_stock_17_04_20.hdf5": {
        "url": "https://ndownloader.figshare.com/files/23086469",
        "size": 663_232_280,
        "md5": "00de71724ec1d4a8463c1b0d8a5d0941",
    },
}
OUT_DIR = Path("data/aizynth")


def build_config_yaml(out_dir: Path) -> str:
    """Render AiZynthFinder `config.yml`, including `stock:` only when the HDF5 stock file exists."""
    resolved = out_dir.resolve()
    lines = [
        "expansion:",
        "  uspto:",
        f"    - {resolved / 'uspto_model.onnx'}",
        f"    - {resolved / 'uspto_templates.csv.gz'}",
        "filter:",
        f"  uspto: {resolved / 'uspto_filter_model.onnx'}",
    ]
    stock_file = out_dir / "zinc_stock_17_04_20.hdf5"
    if stock_file.exists() and stock_file.stat().st_size > 0:
        lines.extend(
            [
                "stock:",
                f"  zinc: {resolved / 'zinc_stock_17_04_20.hdf5'}",
            ]
        )
    return "\n".join(lines) + "\n"


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Fetch AiZynthFinder public data (~790 MB total).")
    parser.add_argument("--with-stock", action="store_true", help="Also download 663 MB ZINC HDF5 stock.")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR, help="Output directory.")
    parser.add_argument("--force", action="store_true", help="Re-download even if verified.")
    args = parser.parse_args(argv)

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    targets = ["uspto_model.onnx", "uspto_templates.csv.gz", "uspto_filter_model.onnx"]
    if args.with_stock:
        targets.append("zinc_stock_17_04_20.hdf5")

    print(
        "Note: If Zenodo returns 403 from automated clients, run `download_public_data data/aizynth` "
        "after `pip install -e .[aizynthfinder]`."
    )

    for name in targets:
        spec = FILES[name]
        dest = out_dir / name
        url = str(spec["url"])
        print(f"Downloading/verifying {url} -> {dest}")
        try:
            sha = atomic_download(
                url,
                dest,
                expected_md5=str(spec["md5"]),
                expected_size=int(spec["size"]),  # type: ignore[arg-type]
                timeout=180.0,
                force=args.force,
            )
            print(f"Verified {dest.name} (sha256={sha})")
        except (httpx.HTTPError, ValueError) as exc:
            print(f"Failed to download/verify {name}: {exc}")
            return 1

    cfg = out_dir / "config.yml"
    cfg.write_text(build_config_yaml(out_dir), encoding="utf-8")
    print(f"Wrote {cfg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
