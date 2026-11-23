#!/usr/bin/env python3
"""Fetch SCScore MIT-licensed NumPy weights (6.18 MB) into data/scscore/ with SHA-256 verification."""

from __future__ import annotations

import argparse
from pathlib import Path

from buyormake.canonical import atomic_download

URL = (
    "https://raw.githubusercontent.com/connorcoley/scscore/master/"
    "models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz"
)
EXPECTED_SHA256 = "fa90bd552307b2efd05eb1e53033bc9eadadfd5bdca682e541040982af99536b"
EXPECTED_SIZE = 6_178_760
OUT_PATH = Path("data/scscore/model.ckpt-10654.as_numpy.json.gz")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Fetch MIT-licensed SCScore NumPy weights.")
    parser.add_argument("--out", type=Path, default=OUT_PATH, help="Destination path.")
    parser.add_argument("--force", action="store_true", help="Re-download even if verified.")
    args = parser.parse_args(argv)

    print(f"Ensuring verified SCScore weights -> {args.out}")
    digest = atomic_download(
        URL,
        args.out,
        expected_sha256=EXPECTED_SHA256,
        expected_size=EXPECTED_SIZE,
        timeout=60.0,
        force=args.force,
    )
    print(f"Verified {args.out} ({args.out.stat().st_size} bytes, sha256={digest})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
