"""Tests for LocalStockBackend multi-format loaders (.csv, .jsonl, .smi.gz, .sqlite, .hdf5)."""

from __future__ import annotations

import gzip
import json
import sqlite3
import sys
import types
from pathlib import Path

import numpy as np
import pytest

from buyormake import cost_aware_rank
from buyormake.canonical import normalize_smiles
from buyormake.stock import LocalStockBackend


def test_csv_jsonl_and_smi_gz_stock_loaders(tmp_path: Path) -> None:
    # 1. Priced .csv
    csv_path = tmp_path / "catalog.csv"
    csv_path.write_text(
        "smiles,vendor,catalog_id,pack_value,pack_unit,price,currency,lead_time_days,in_stock,source_id,source_url,license\n"
        "CCO,CsvVendor,CSV-1,10,g,25.0,USD,3,true,local_stock,https://example.org/csv1,CC0-1.0\n",
        encoding="utf-8",
    )

    # 2. Priced .jsonl
    prop_nm = normalize_smiles("CCCO")
    assert prop_nm is not None
    jsonl_path = tmp_path / "catalog.jsonl"
    jsonl_path.write_text(
        json.dumps(
            {
                "inchikey": prop_nm.inchikey,
                "vendor": "JsonlVendor",
                "catalog_id": "JL-1",
                "pack_size": {"value": 500.0, "unit": "mg", "raw_unit_str": "500 mg"},
                "price": 10.0,
                "currency": "USD",
                "price_basis": "listed_catalog",
                "lead_time_days": 5,
                "in_stock": True,
                "source_id": "local_stock",
                "source_url": "https://example.org/jl1",
                "license": "CC0-1.0",
            }
        )
        + "\n",
        encoding="utf-8",
    )

    # 3. Compressed .smi.gz
    smigz_path = tmp_path / "building_blocks.smi.gz"
    with gzip.open(smigz_path, "wt", encoding="utf-8") as gz:
        gz.write("# compressed SMILES\nc1ccccc1 BENZ-01\n")

    ranked = cost_aware_rank(
        ["CCO", "CCCO", "c1ccccc1"],
        max_steps=1,
        backends=("local_stock",),
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[csv_path, jsonl_path, smigz_path],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    # CCO = $25/10g = 2.5 USD/g; CCCO = $10/0.5g = 20.0 USD/g; benzene = unpriced
    assert [r.candidate for r in ranked] == ["CCO", "CCCO", "c1ccccc1"]
    assert ranked[0].best_price is not None and ranked[0].best_price.value == pytest.approx(2.5)
    assert ranked[1].best_price is not None and ranked[1].best_price.value == pytest.approx(20.0)
    assert ranked[2].best_price is None


def test_sqlite_priced_and_legacy_stock_columns(tmp_path: Path) -> None:
    nm_priced = normalize_smiles("CC(=O)Oc1ccccc1C(=O)O")
    nm_legacy = normalize_smiles("CCO")
    assert nm_priced and nm_legacy

    # Priced SQLite catalog
    priced_db = tmp_path / "priced_stock.sqlite"
    conn = sqlite3.connect(str(priced_db))
    conn.execute(
        """
        CREATE TABLE stock (
            inchikey TEXT NOT NULL,
            vendor TEXT NOT NULL,
            catalog_id TEXT,
            pack_value REAL,
            pack_unit TEXT,
            purity TEXT,
            price REAL,
            currency TEXT,
            price_basis TEXT,
            lead_time_days INTEGER,
            in_stock INTEGER,
            source_id TEXT NOT NULL,
            source_url TEXT NOT NULL,
            fetched_utc TEXT NOT NULL,
            license TEXT NOT NULL
        )
        """
    )
    conn.execute(
        "INSERT INTO stock VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (
            nm_priced.inchikey,
            "SqliteChem",
            "SQ-2244",
            25.0,
            "g",
            "99%",
            50.0,
            "USD",
            "listed_catalog",
            2,
            1,
            "local_stock",
            "https://example.org/sq2244",
            "2026-11-18T00:00:00Z",
            "CC0-1.0",
        ),
    )
    conn.commit()
    conn.close()

    # Legacy 7-column unpriced SQLite catalog
    legacy_db = tmp_path / "legacy_stock.sqlite"
    conn2 = sqlite3.connect(str(legacy_db))
    conn2.execute(
        "CREATE TABLE stock (inchikey TEXT, vendor TEXT, catalog_id TEXT, "
        "source_id TEXT, source_url TEXT, fetched_utc TEXT, license TEXT)"
    )
    conn2.execute(
        "INSERT INTO stock VALUES (?, ?, ?, ?, ?, ?, ?)",
        (
            nm_legacy.inchikey,
            "LegacyVendor",
            "LG-1",
            "local_stock",
            "https://example.org/lg1",
            "2026-11-18T00:00:00Z",
            "CC0-1.0",
        ),
    )
    conn2.commit()
    conn2.close()

    backend = LocalStockBackend([priced_db, legacy_db])
    try:
        offs = backend.lookup(nm_priced.inchikey)
        assert len(offs) == 1
        assert offs[0].price == 50.0
        assert offs[0].currency == "USD"
        assert offs[0].price_basis == "listed_catalog"
        assert offs[0].pack_size is not None and offs[0].pack_size.value == 25.0
        assert offs[0].lead_time_days == 2
        assert offs[0].in_stock is True

        offs_leg = backend.lookup(nm_legacy.inchikey)
        assert len(offs_leg) == 1
        assert offs_leg[0].price is None
        assert offs_leg[0].price_basis == "unknown"
    finally:
        backend.close()

    ranked = cost_aware_rank(
        ["CC(=O)Oc1ccccc1C(=O)O"],
        backends=("local_stock",),
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[priced_db],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert ranked[0].best_price is not None
    assert ranked[0].best_price.unit == "USD/g"
    assert ranked[0].best_price.value == pytest.approx(2.0)


def test_hdf5_stock_warns_and_raises_on_unreadable_and_supports_h5py_fallback(
    tmp_path: Path, monkeypatch
) -> None:
    bad_h5 = tmp_path / "corrupt_or_missing_deps.hdf5"
    bad_h5.write_bytes(b"not a valid hdf5 file")

    with pytest.warns(UserWarning, match="Failed to read HDF5 stock file"), pytest.raises(
        RuntimeError, match="Failed to read HDF5 stock file"
    ):
        LocalStockBackend([bad_h5])

    # Verify h5py fallback reader on PyTables-style structured table dataset
    nm = normalize_smiles("CCO")
    assert nm is not None

    class FakeDataset:
        def __init__(self, arr: np.ndarray) -> None:
            self._arr = arr
            self.dtype = arr.dtype
            self.shape = arr.shape

        def __getitem__(self, key):
            return self._arr[key]

    struct_arr = np.array(
        [(nm.inchikey.encode("utf-8"),)],
        dtype=[("inchi_key", "S27")],
    )

    class FakeH5File:
        def __init__(self, path, mode="r") -> None:
            self._groups = {"table": {"table": FakeDataset(struct_arr)}}

        def __enter__(self):
            return self._groups

        def __exit__(self, *args):
            return None

    fake_h5py = types.ModuleType("h5py")
    fake_h5py.File = FakeH5File  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "h5py", fake_h5py)

    good_h5 = tmp_path / "zinc_mock.hdf5"
    good_h5.write_bytes(b"mock")
    backend = LocalStockBackend([good_h5])
    found = backend.lookup(nm.inchikey)
    assert len(found) == 1
    assert found[0].source_id == "aizynthfinder_zinc_stock"
