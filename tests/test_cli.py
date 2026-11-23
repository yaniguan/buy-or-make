"""Tests for the buy-or-make CLI."""

from __future__ import annotations

import json
from pathlib import Path

from buyormake.cli import main
from buyormake.schema import AvailabilityRecord


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


def test_cli_names_records_out_and_budget_filter(tmp_path: Path, capsys) -> None:
    smiles_file = tmp_path / "candidates.smi"
    smiles_file.write_text(
        "# SMILES Name\n"
        "CC(=O)Nc1ccc(O)cc1 Paracetamol\n"
        "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
        "CCO EthanolUnpriced\n",
        encoding="utf-8",
    )
    stock_tsv = tmp_path / "stock.tsv"
    stock_tsv.write_text(
        "smiles\tvendor\tcatalog_id\tpack_value\tpack_unit\tprice\tcurrency\tsource_id\tsource_url\tlicense\n"
        "CC(=O)Nc1ccc(O)cc1\tVendorA\tP-1\t5\tg\t10.0\tUSD\tlocal_stock\thttps://example.org/p1\tCC0-1.0\n"
        "CC(=O)Oc1ccccc1C(=O)O\tVendorA\tA-1\t1\tg\t50.0\tUSD\tlocal_stock\thttps://example.org/a1\tCC0-1.0\n"
        "CCO\tVendorB\tE-1\t\t\t\t\tlocal_stock\thttps://example.org/e1\tCC0-1.0\n",
        encoding="utf-8",
    )
    records_jsonl = tmp_path / "records.jsonl"

    rc = main(
        [
            "filter",
            str(smiles_file),
            "--backends",
            "local_stock",
            "--offline",
            "--stock",
            str(stock_tsv),
            "--records-out",
            str(records_jsonl),
            "--budget",
            "10.0",
            "--require-complete-cost",
            "--json",
        ]
    )
    assert rc == 0
    summary = json.loads(capsys.readouterr().out)
    assert [r["name"] for r in summary["kept"]] == ["Paracetamol"]
    assert [r["name"] for r in summary["over_budget"]] == ["Aspirin"]
    assert len(summary["rejected"]) == 2  # Aspirin (over_budget) + EthanolUnpriced (missing price)

    lines = [ln for ln in records_jsonl.read_text(encoding="utf-8").splitlines() if ln.strip()]
    assert len(lines) == 3
    recs = [AvailabilityRecord.model_validate_json(ln) for ln in lines]
    assert [r.molecule.name for r in recs] == ["Paracetamol", "Aspirin", "EthanolUnpriced"]
    assert all(r.record_hash == r.content_hash() for r in recs)
