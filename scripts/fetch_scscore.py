#!/usr/bin/env python3
"""Fetch SCScore MIT-licensed NumPy weights (6.18 MB) into data/scscore/."""

from __future__ import annotations

from pathlib import Path
import httpx

URL = (
    "https://raw.githubusercontent.com/connorcoley/scscore/master/"
    "models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz"
)
OUT_PATH = Path("data/scscore/model.ckpt-10654.as_numpy.json.gz")


def main() -> int:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    if OUT_PATH.exists() and OUT_PATH.stat().st_size > 1_000_000:
        print(f"Already present: {OUT_PATH} ({OUT_PATH.stat().st_size} bytes)")
        return 0
    print(f"Downloading {URL} -> {OUT_PATH}")
    with httpx.stream("GET", URL, follow_redirects=True, timeout=60.0) as r:
        r.raise_for_status()
        with OUT_PATH.open("wb") as f:
            for chunk in r.iter_bytes():
                f.write(chunk)
    print(f"Saved {OUT_PATH} ({OUT_PATH.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
