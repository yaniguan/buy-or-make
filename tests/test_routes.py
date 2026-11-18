"""Tests for retrosynthesis routes: 1-step/2-step template fallback and SQLite reaction lookup."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from buyormake import cost_aware_rank
from buyormake.canonical import normalize_smiles


def test_reaction_lookup_backend(tmp_path: Path) -> None:
    db_path = tmp_path / "reactions.sqlite"
    conn = sqlite3.connect(str(db_path))
    conn.execute(
        "CREATE TABLE reactions (id INTEGER PRIMARY KEY, product_inchikey TEXT, "
        "product_smiles TEXT, reactants_json TEXT, reaction_smiles TEXT, reference TEXT, source TEXT)"
    )
    target = "CC(=O)Nc1ccc(F)cc1"
    nm = normalize_smiles(target)
    assert nm is not None
    conn.execute(
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, reaction_smiles, reference, source) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        (
            nm.inchikey,
            nm.canonical_smiles,
            json.dumps(["CC(=O)O", "Nc1ccc(F)cc1"]),
            "CC(=O)O.Nc1ccc(F)cc1>>CC(=O)Nc1ccc(F)cc1",
            "US1234567",
            "uspto_lookup",
        ),
    )
    conn.commit()
    conn.close()

    stock = tmp_path / "stock.smi"
    stock.write_text("CC(=O)O acetic_acid\nNc1ccc(F)cc1 4_fluoroaniline\n", encoding="utf-8")

    res = cost_aware_rank(
        [target],
        max_steps=1,
        backends=("local_stock",),
        route_backend="lookup",
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock],
        reaction_db_path=db_path,
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert len(res) == 1
    assert res[0].status == "makeable"
    assert res[0].n_steps == 1
    assert any("uspto_lookup" in reason for reason in res[0].reasons)
