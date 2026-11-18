#!/usr/bin/env python3
"""Fetch Open Reaction Database (ord-data, CC-BY-SA-4.0) parquet shards and index into data/reactions.sqlite."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

import httpx

from buyormake.canonical import normalize_smiles

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
            reaction_smiles TEXT NOT NULL,
            reference TEXT NOT NULL,
            source TEXT NOT NULL
        )
        """
    )
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rxn_prod_ik ON reactions(product_inchikey)")
    conn.commit()


def list_shards(shard_prefix: str = "00") -> list[str]:
    r = httpx.get(f"{HF_API_TREE}/{shard_prefix}", timeout=20.0)
    r.raise_for_status()
    return [item["path"] for item in r.json() if item.get("path", "").endswith(".parquet")]


def download_shard(repo_path: str, out_dir: Path = OUT_DIR) -> Path:
    dest = out_dir / Path(repo_path).name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = f"{HF_RESOLVE}/{repo_path}"
    print(f"Downloading {url} -> {dest}")
    with httpx.stream("GET", url, follow_redirects=True, timeout=120.0) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_bytes():
                f.write(chunk)
    return dest


def index_parquet(parquet_path: Path, db_path: Path = DB_PATH, max_rows: int = 50_000) -> int:
    try:
        import pyarrow.parquet as pq  # type: ignore
    except ImportError:
        print("Install pyarrow (pip install -e .[data]) to parse ORD parquet files.")
        return 0

    table = pq.read_table(parquet_path)
    cols = set(table.column_names)
    smi_col = next((c for c in ("reaction_smiles", "rxn_smiles", "smiles") if c in cols), None)
    id_col = next((c for c in ("reaction_id", "ord_id", "dataset_id") if c in cols), None)
    if smi_col is None:
        print(f"No reaction SMILES column found in {parquet_path.name}; columns={sorted(cols)}")
        return 0

    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    _init_db(conn)
    inserted = 0
    batch = []
    smis = table[smi_col].to_pylist()
    ids = table[id_col].to_pylist() if id_col else [f"ord:{i}" for i in range(len(smis))]
    for rxn_smi, ref in zip(smis, ids):
        if not rxn_smi or ">>" not in rxn_smi:
            continue
        lhs, rhs = rxn_smi.split(">>", 1)
        reactants = [r.strip() for r in lhs.split(">")[0].split(".") if r.strip()]
        products = [p.strip() for p in rhs.split(".") if p.strip()]
        for p_raw in products:
            nm = normalize_smiles(p_raw)
            if nm is None:
                continue
            batch.append(
                (
                    nm.inchikey,
                    nm.canonical_smiles,
                    json.dumps(reactants),
                    rxn_smi,
                    str(ref or "ord-data"),
                    "ord_lookup",
                )
            )
            inserted += 1
            if max_rows and inserted >= max_rows:
                break
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
    parser = argparse.ArgumentParser(description="Fetch ORD ord-data parquet shards (CC-BY-SA-4.0).")
    parser.add_argument("--prefix", default="00", help="Shard directory prefix under data/ (default: 00).")
    parser.add_argument("--max-files", type=int, default=1, help="Number of parquet files to download.")
    parser.add_argument("--max-rows", type=int, default=50_000, help="Max rows per file to index.")
    args = parser.parse_args()

    shards = list_shards(args.prefix)[: args.max_files]
    for shard in shards:
        p = download_shard(shard)
        n = index_parquet(p, DB_PATH, max_rows=args.max_rows)
        print(f"Indexed {n} rows from {p.name} -> {DB_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
