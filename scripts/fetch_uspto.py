#!/usr/bin/env python3
"""Fetch Lowe's USPTO reactions (CC0-1.0, Figshare 5104873) and build data/reactions.sqlite."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Optional

from buyormake.canonical import atomic_download
from buyormake.routes import parse_reaction_smiles_for_indexing

USPTO_GRANTS_URL = "https://ndownloader.figshare.com/files/8664379"  # 1976_Sep2016_USPTOgrants_smiles.7z
USPTO_GRANTS_SIZE = 75_191_394
USPTO_GRANTS_MD5 = "b0d94806fe7479fad6db9b894cff858e"
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
            agents_json TEXT NOT NULL DEFAULT '[]',
            reaction_smiles TEXT NOT NULL,
            reference TEXT NOT NULL,
            source TEXT NOT NULL
        )
        """
    )
    cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(reactions)").fetchall()}
    if "agents_json" not in cols:
        conn.execute("ALTER TABLE reactions ADD COLUMN agents_json TEXT NOT NULL DEFAULT '[]'")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rxn_prod_ik ON reactions(product_inchikey)")
    conn.commit()


def _download(
    url: str,
    dest: Path,
    *,
    expected_sha256: Optional[str] = None,
    expected_md5: Optional[str] = USPTO_GRANTS_MD5,
    expected_size: Optional[int] = USPTO_GRANTS_SIZE,
    force: bool = False,
) -> str:
    print(f"Downloading/verifying {url} -> {dest}")
    return atomic_download(
        url,
        dest,
        expected_sha256=expected_sha256,
        expected_md5=expected_md5,
        expected_size=expected_size,
        timeout=180.0,
        force=force,
    )


def index_smi_file(smi_path: Path, db_path: Path = DB_PATH, max_rows: int = 100_000) -> int:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    _init_db(conn)
    inserted = 0
    batch: list[tuple[str, str, str, str, str, str, str]] = []
    with smi_path.open("r", encoding="utf-8", errors="replace") as f:
        for idx, line in enumerate(f):
            if idx == 0 and "ReactionSmiles" in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            rxn_smi = parts[0].strip()
            patent_id = parts[1].strip() if len(parts) > 1 and parts[1].strip() else f"uspto:{idx}"
            for p_nm, stoich_reactants, agents in parse_reaction_smiles_for_indexing(rxn_smi):
                batch.append(
                    (
                        p_nm.inchikey,
                        p_nm.canonical_smiles,
                        json.dumps(stoich_reactants),
                        json.dumps(agents),
                        rxn_smi,
                        patent_id,
                        "uspto_lookup",
                    )
                )
                inserted += 1
                if len(batch) >= 2000:
                    conn.executemany(
                        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
                        "agents_json, reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?, ?)",
                        batch,
                    )
                    conn.commit()
                    batch.clear()
                if max_rows and inserted >= max_rows:
                    break
            if max_rows and inserted >= max_rows:
                break
    if batch:
        conn.executemany(
            "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
            "agents_json, reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?, ?)",
            batch,
        )
        conn.commit()
    conn.close()
    return inserted


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Fetch Lowe USPTO reaction SMILES (CC0-1.0).")
    parser.add_argument("--max-rows", type=int, default=100_000, help="Max product rows to index (0 = all).")
    parser.add_argument("--expected-sha256", default=None, help="Optional pinned SHA-256 hex digest.")
    parser.add_argument("--force", action="store_true", help="Re-download archive even if verified.")
    args = parser.parse_args(argv)

    archive = OUT_DIR / "1976_Sep2016_USPTOgrants_smiles.7z"
    sha = _download(USPTO_GRANTS_URL, archive, expected_sha256=args.expected_sha256, force=args.force)
    print(f"Verified archive sha256={sha}")

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
