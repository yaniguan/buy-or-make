#!/usr/bin/env python3
"""Fetch Lowe's USPTO reactions (CC0-1.0, Figshare 5104873) and build data/reactions.sqlite."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

import httpx

from buyormake.canonical import normalize_smiles

USPTO_GRANTS_URL = "https://ndownloader.figshare.com/files/8664379"  # 1976_Sep2016_USPTOgrants_smiles.7z (75.2 MB)
OUT_DIR = Path("data/uspto")
DB_PATH = Path("data/reactions.sqlite")


def _init_db(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS reactions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            product_inchikey TEXT NOT NULL,
            product_smiles TEXT NOT NULL,
            reactants_json TEXT NOT NULL,
            reaction_smiles TEXT NOT NULL,
            reference TEXT NOT NULL,
            source TEXT NOT NULL
        )
        """
    )
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rxn_prod_ik ON reactions(product_inchikey)")
    conn.commit()


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1_000_000:
        print(f"Already downloaded: {dest}")
        return
    print(f"Downloading {url} -> {dest}")
    with httpx.stream("GET", url, follow_redirects=True, timeout=120.0) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_bytes():
                f.write(chunk)


def index_smi_file(smi_path: Path, db_path: Path = DB_PATH, max_rows: int = 100_000) -> int:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    _init_db(conn)
    inserted = 0
    batch = []
    with smi_path.open("r", encoding="utf-8", errors="replace") as f:
        for idx, line in enumerate(f):
            if idx == 0 and "ReactionSmiles" in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            rxn_smi = parts[0].strip()
            patent_id = parts[1].strip() if len(parts) > 1 else f"uspto:{idx}"
            if ">>" not in rxn_smi:
                continue
            lhs, rhs = rxn_smi.split(">>", 1)
            reactant_part = lhs.split(">")[0]
            reactants = [r.strip() for r in reactant_part.split(".") if r.strip()]
            products = [p.strip() for p in rhs.split(".") if p.strip()]
            if not reactants or not products:
                continue
            for p_raw in products:
                p_nm = normalize_smiles(p_raw.split()[0])
                if p_nm is None:
                    continue
                batch.append(
                    (
                        p_nm.inchikey,
                        p_nm.canonical_smiles,
                        json.dumps(reactants),
                        rxn_smi,
                        patent_id,
                        "uspto_lookup",
                    )
                )
                inserted += 1
                if len(batch) >= 2000:
                    conn.executemany(
                        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
                        "reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?)",
                        batch,
                    )
                    conn.commit()
                    batch.clear()
            if max_rows and inserted >= max_rows:
                break
    if batch:
        conn.executemany(
            "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
            "reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?)",
            batch,
        )
        conn.commit()
    conn.close()
    return inserted


def main() -> int:
    parser = argparse.ArgumentParser(description="Fetch Lowe USPTO reaction SMILES (CC0-1.0).")
    parser.add_argument("--max-rows", type=int, default=100_000, help="Max product rows to index (0 = all).")
    args = parser.parse_args()

    archive = OUT_DIR / "1976_Sep2016_USPTOgrants_smiles.7z"
    _download(USPTO_GRANTS_URL, archive)

    try:
        import py7zr  # type: ignore
    except ImportError:
        print("Install py7zr (pip install -e .[data]) to extract .7z archives.")
        return 1

    with py7zr.SevenZipFile(archive, mode="r") as z:
        z.extractall(path=OUT_DIR)

    smi_files = list(OUT_DIR.glob("*.rsmi")) + list(OUT_DIR.glob("*.smi")) + list(OUT_DIR.glob("*.txt"))
    if not smi_files:
        print(f"No extracted reaction files found in {OUT_DIR}")
        return 1

    total = 0
    for f in smi_files:
        n = index_smi_file(f, DB_PATH, max_rows=args.max_rows)
        print(f"Indexed {n} rows from {f.name} -> {DB_PATH}")
        total += n
    print(f"Done. Total indexed: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
