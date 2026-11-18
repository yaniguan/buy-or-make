"""Tests for the buy-or-make CLI."""

from __future__ import annotations

import json
from pathlib import Path

from buyormake.cli import main


def test_cli_rank_json(tmp_path: Path, capsys) -> None:
    smiles_file = tmp_path / "smiles.txt"
    smiles_file.write_text("CCO\nbad_smiles\n", encoding="utf-8")
    stock_file = tmp_path / "stock.smi"
    stock_file.write_text("CCO ethanol\n", encoding="utf-8")

    rc = main(
        [
            "rank",
            str(smiles_file),
            "--max-steps",
            "2",
            "--backends",
            "local_stock",
            "--offline",
            "--stock",
            str(stock_file),
            "--json",
        ]
    )
    assert rc == 0
    out = capsys.readouterr().out
    data = json.loads(out)
    assert [row["status"] for row in data] == ["buyable", "invalid"]
