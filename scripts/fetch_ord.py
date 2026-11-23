#!/usr/bin/env python3
"""Fetch Open Reaction Database (ord-data, CC-BY-SA-4.0) parquet shards and index into data/reactions.sqlite."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Iterable, Optional

import httpx

from buyormake.canonical import atomic_download
from buyormake.routes import parse_reaction_smiles_for_indexing

HF_API_TREE = "https://huggingface.co/api/datasets/open-reaction-database/ord-data/tree/main/data"
HF_RESOLVE = "https://huggingface.co/datasets/open-reaction-database/ord-data/resolve/main"
OUT_DIR = Path("data/ord")
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


def list_shards_with_checksums(shard_prefix: str = "00") -> list[tuple[str, Optional[str], Optional[int]]]:
    """Return `(repo_path, lfs_sha256, byte_size)` from Hugging Face Hub metadata."""
    r = httpx.get(f"{HF_API_TREE}/{shard_prefix}", timeout=20.0)
    r.raise_for_status()
    out: list[tuple[str, Optional[str], Optional[int]]] = []
    for item in r.json():
        path = item.get("path", "")
        if not path.endswith(".parquet"):
            continue
        lfs = item.get("lfs") or {}
        sha256 = lfs.get("oid") if isinstance(lfs.get("oid"), str) and len(lfs.get("oid")) == 64 else None
        size = lfs.get("size") or item.get("size")
        size_int = int(size) if isinstance(size, int) and size > 0 else None
        out.append((path, sha256, size_int))
    return out


def list_shards(shard_prefix: str = "00") -> list[str]:
    return [path for path, _, _ in list_shards_with_checksums(shard_prefix)]


def download_shard(
    repo_path: str,
    out_dir: Path = OUT_DIR,
    *,
    expected_sha256: Optional[str] = None,
    expected_size: Optional[int] = None,
    force: bool = False,
) -> Path:
    dest = out_dir / Path(repo_path).name
    url = f"{HF_RESOLVE}/{repo_path}"
    print(f"Downloading/verifying {url} -> {dest}")
    sha = atomic_download(
        url,
        dest,
        expected_sha256=expected_sha256,
        expected_size=expected_size,
        timeout=120.0,
        force=force,
    )
    print(f"Verified {dest.name} (sha256={sha})")
    return dest


def index_ord_rows(
    rows: Iterable[tuple[Optional[str], Optional[str]]],
    db_path: Path = DB_PATH,
    max_rows: int = 50_000,
) -> int:
    """Index `(reaction_smiles, reaction_id)` pairs into `reactions.sqlite`, filtering spectators."""
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    _init_db(conn)
    inserted = 0
    batch: list[tuple[str, str, str, str, str, str, str]] = []
    for idx, (rxn_smi, ref) in enumerate(rows):
        if not rxn_smi or ">>" not in rxn_smi:
            continue
        ref_str = str(ref or f"ord:{idx}")
        for p_nm, stoich_reactants, agents in parse_reaction_smiles_for_indexing(rxn_smi):
            batch.append(
                (
                    p_nm.inchikey,
                    p_nm.canonical_smiles,
                    json.dumps(stoich_reactants),
                    json.dumps(agents),
                    rxn_smi,
                    ref_str,
                    "ord_lookup",
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


def index_parquet(parquet_path: Path, db_path: Path = DB_PATH, max_rows: int = 50_000) -> int:
    smis: list[Optional[str]]
    ids: list[Optional[str]]
    try:
        import pyarrow.parquet as pq  # type: ignore

        table = pq.read_table(parquet_path)
        cols = set(table.column_names)
        smi_col = next((c for c in ("reaction_smiles", "rxn_smiles", "smiles") if c in cols), None)
        id_col = next((c for c in ("reaction_id", "ord_id", "dataset_id") if c in cols), None)
        if smi_col is None:
            print(f"No reaction SMILES column found in {parquet_path.name}; columns={sorted(cols)}")
            return 0
        smis = table[smi_col].to_pylist()
        ids = table[id_col].to_pylist() if id_col else [f"ord:{i}" for i in range(len(smis))]
    except ImportError:
        try:
            import pandas as pd  # type: ignore

            df = pd.read_parquet(parquet_path)
            cols = set(df.columns)
            smi_col = next((c for c in ("reaction_smiles", "rxn_smiles", "smiles") if c in cols), None)
            id_col = next((c for c in ("reaction_id", "ord_id", "dataset_id") if c in cols), None)
            if smi_col is None:
                print(f"No reaction SMILES column found in {parquet_path.name}; columns={sorted(cols)}")
                return 0
            smis = df[smi_col].tolist()
            ids = df[id_col].tolist() if id_col else [f"ord:{i}" for i in range(len(smis))]
        except Exception:
            print("Install pyarrow (pip install -e .[data]) to parse ORD parquet files.")
            return 0

    return index_ord_rows(zip(smis, ids), db_path=db_path, max_rows=max_rows)


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Fetch ORD ord-data parquet shards (CC-BY-SA-4.0).")
    parser.add_argument("--prefix", default="00", help="Shard directory prefix under data/ (default: 00).")
    parser.add_argument("--max-files", type=int, default=1, help="Number of parquet files to download.")
    parser.add_argument("--max-rows", type=int, default=50_000, help="Max rows per file to index.")
    parser.add_argument("--force", action="store_true", help="Re-download shards even if verified.")
    args = parser.parse_args(argv)

    shards = list_shards_with_checksums(args.prefix)[: args.max_files]
    for repo_path, sha256, size in shards:
        p = download_shard(
            repo_path,
            OUT_DIR,
            expected_sha256=sha256,
            expected_size=size,
            force=args.force,
        )
        n = index_parquet(p, DB_PATH, max_rows=args.max_rows)
        print(f"Indexed {n} rows from {p.name} -> {DB_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
