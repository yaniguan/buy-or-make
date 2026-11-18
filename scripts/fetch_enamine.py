#!/usr/bin/env python3
"""Manual guide & local SDF importer for Enamine Building Blocks (requires registration; fetch_only)."""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path
from rdkit import Chem

from buyormake.canonical import normalize_smiles
from buyormake.schema import utc_now

DB_PATH = Path("data/stock.sqlite")


def import_sdf(sdf_path: Path, db_path: Path = DB_PATH, max_rows: int = 0) -> int:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS stock (
            inchikey TEXT NOT NULL,
            vendor TEXT NOT NULL,
            catalog_id TEXT,
            source_id TEXT NOT NULL,
            source_url TEXT NOT NULL,
            fetched_utc TEXT NOT NULL,
            license TEXT NOT NULL
        )
        """
    )
    conn.execute("CREATE INDEX IF NOT EXISTS idx_stock_ik ON stock(inchikey)")
    ts = utc_now()
    suppl = Chem.SDMolSupplier(str(sdf_path))
    batch = []
    count = 0
    for mol in suppl:
        if mol is None:
            continue
        smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        nm = normalize_smiles(smi)
        if nm is None:
            continue
        cat_id = None
        for prop in ("Catalog_ID", "ID", "ENAMINE_ID", "_Name"):
            if prop == "_Name" and mol.HasProp("_Name"):
                cat_id = mol.GetProp("_Name").strip() or None
                break
            if mol.HasProp(prop):
                cat_id = mol.GetProp(prop).strip() or None
                break
        batch.append(
            (
                nm.inchikey,
                "Enamine",
                cat_id,
                "enamine_bb",
                "https://enamine.net/building-blocks/building-blocks-catalog",
                ts,
                "custom/terms",
            )
        )
        count += 1
        if len(batch) >= 2000:
            conn.executemany("INSERT INTO stock VALUES (?, ?, ?, ?, ?, ?, ?)", batch)
            conn.commit()
            batch.clear()
        if max_rows and count >= max_rows:
            break
    if batch:
        conn.executemany("INSERT INTO stock VALUES (?, ?, ?, ?, ?, ?, ?)", batch)
        conn.commit()
    conn.close()
    return count


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Import a manually downloaded Enamine Building Blocks SDF into data/stock.sqlite."
    )
    parser.add_argument("--sdf", type=Path, default=None, help="Path to locally downloaded Enamine SDF.")
    parser.add_argument("--max-rows", type=int, default=0, help="Optional limit for testing.")
    args = parser.parse_args()

    if args.sdf is None:
        print(
            "Enamine Building Blocks require user registration and acceptance of Enamine Terms of Use:\n"
            "  1. Visit https://enamine.net/building-blocks/building-blocks-catalog\n"
            "  2. Register / sign in and download the desired building-block SDF to data/enamine/\n"
            "  3. Run: python scripts/fetch_enamine.py --sdf data/enamine/<file>.sdf\n"
            "Note: Open Structure Data includes structures + Enamine IDs only (no prices), and "
            "Enamine Terms of Use forbid republication or commercial dataset redistribution."
        )
        return 0

    if not args.sdf.exists():
        print(f"File not found: {args.sdf}")
        return 1
    n = import_sdf(args.sdf, max_rows=args.max_rows)
    print(f"Imported {n} Enamine structures into {DB_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
