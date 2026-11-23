"""Tests for reaction/SDF indexers, atomic checksum downloads, and AiZynthFinder config generation."""

from __future__ import annotations

import hashlib
import importlib.util
import sqlite3
from pathlib import Path

import httpx
import pytest
from rdkit import Chem

from buyormake.canonical import atomic_download, normalize_smiles
from buyormake.routes import ReactionLookupBackend, _select_if_configured
from buyormake.stock import LocalStockBackend

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_script_module(script_name: str):
    path = REPO_ROOT / "scripts" / script_name
    spec = importlib.util.spec_from_file_location(script_name[:-3], path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_enamine_import_sdf_into_sqlite(tmp_path: Path) -> None:
    fetch_enamine = _load_script_module("fetch_enamine.py")
    sdf_path = tmp_path / "enamine_sample.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    m1 = Chem.MolFromSmiles("Nc1ccc(Br)cc1")
    m2 = Chem.MolFromSmiles("OB(O)c1ccccn1")
    assert m1 and m2
    m1.SetProp("ENAMINE_ID", "EN300-001")
    m2.SetProp("_Name", "EN300-002")
    writer.write(m1)
    writer.write(m2)
    writer.close()

    db_path = tmp_path / "stock.sqlite"
    count = fetch_enamine.import_sdf(sdf_path, db_path=db_path)
    assert count == 2

    backend = LocalStockBackend([db_path])
    try:
        nm1 = normalize_smiles("Nc1ccc(Br)cc1")
        assert nm1 is not None
        offs = backend.lookup(nm1.inchikey)
        assert len(offs) == 1
        assert offs[0].vendor == "Enamine"
        assert offs[0].catalog_id == "EN300-001"
        assert offs[0].source_id == "enamine_bb"
    finally:
        backend.close()


def test_uspto_and_ord_indexers_strip_spectators(tmp_path: Path) -> None:
    fetch_uspto = _load_script_module("fetch_uspto.py")
    fetch_ord = _load_script_module("fetch_ord.py")

    rsmi_path = tmp_path / "sample.rsmi"
    # Header + mapped Suzuki reaction containing unmapped THF + Pd(PPh3)4
    rxn_line = (
        "[Br:99][c:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1."
        "[OH:97][B:98]([OH:96])[c:7]1[n:8][cH:9][cH:10][cH:11][cH:12]1."
        "C1CCOC1.[Pd]\tUS2020123456A1\n"
    )
    prod = ">>[c:1]1([c:7]2[n:8][cH:9][cH:10][cH:11][cH:12]2)[cH:2][cH:3][cH:4][cH:5][cH:6]1"
    rsmi_path.write_text(
        "ReactionSmiles\tPatentNumber\n" + rxn_line.replace("\t", prod + "\t"),
        encoding="utf-8",
    )

    rxn_db = tmp_path / "reactions.sqlite"
    n_uspto = fetch_uspto.index_smi_file(rsmi_path, db_path=rxn_db)
    assert n_uspto == 1

    # Also index an ORD amide reaction containing DMF solvent in reactant field
    ord_rxn = (
        "[CH3:1][C:2](=[O:3])[OH:90].[NH2:4][c:5]1[cH:6][cH:7][c:8]([F:9])[cH:10][cH:11]1.CN(C)C=O"
        ">>[CH3:1][C:2](=[O:3])[NH:4][c:5]1[cH:6][cH:7][c:8]([F:9])[cH:10][cH:11]1"
    )
    n_ord = fetch_ord.index_ord_rows([(ord_rxn, "ord_dataset-0001")], db_path=rxn_db)
    assert n_ord == 1

    lookup = ReactionLookupBackend(rxn_db)
    try:
        p1 = normalize_smiles("c1ccc(-c2ccccn2)cc1")
        p2 = normalize_smiles("CC(=O)Nc1ccc(F)cc1")
        assert p1 and p2
        st1 = lookup.one_step(p1.inchikey, p1.canonical_smiles)
        assert len(st1) == 1
        assert set(st1[0].reactants) == {"Brc1ccccc1", "OB(O)c1ccccn1"}

        st2 = lookup.one_step(p2.inchikey, p2.canonical_smiles)
        assert len(st2) == 1
        assert set(st2[0].reactants) == {"CC(=O)O", "Nc1ccc(F)cc1"}
    finally:
        lookup.close()

    # Verify agents_json persisted THF and DMF spectators
    conn = sqlite3.connect(str(rxn_db))
    agents_rows = [r[0] for r in conn.execute("SELECT agents_json FROM reactions").fetchall()]
    conn.close()
    assert any("C1CCOC1" in a for a in agents_rows)
    assert any("CN(C)C=O" in a for a in agents_rows)


def test_fetch_aizynth_config_stanza_and_empty_stock_guard(tmp_path: Path) -> None:
    fetch_aizynth = _load_script_module("fetch_aizynth.py")

    # Without --with-stock: zinc_stock_17_04_20.hdf5 is absent -> no `stock:` stanza!
    cfg_no_stock = fetch_aizynth.build_config_yaml(tmp_path)
    assert "expansion:" in cfg_no_stock
    assert "filter:" in cfg_no_stock
    assert "stock:" not in cfg_no_stock

    # With stock file present: `stock:` stanza is included
    (tmp_path / "zinc_stock_17_04_20.hdf5").write_bytes(b"hdf5_placeholder")
    cfg_with_stock = fetch_aizynth.build_config_yaml(tmp_path)
    assert "stock:" in cfg_with_stock
    assert "zinc_stock_17_04_20.hdf5" in cfg_with_stock

    # Guard `_select_if_configured` on empty stock container
    class EmptyStock:
        items: tuple[str, ...] = ()

        def select_all(self) -> None:
            raise AssertionError("select_all() should not be called when items is empty")

    _select_if_configured(EmptyStock())


def test_atomic_download_verifies_sha256_replaces_corrupt_and_cleans_part(tmp_path: Path) -> None:
    payload = b"verified model weights payload" * 100
    good_sha = hashlib.sha256(payload).hexdigest()
    dest = tmp_path / "weights.bin"

    # Pre-seed a corrupted partial file (simulating an interrupted download)
    dest.write_bytes(b"truncated_corrupted_prefix")

    def good_handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(
            200,
            content=payload,
            headers={"Content-Length": str(len(payload))},
        )

    client = httpx.Client(transport=httpx.MockTransport(good_handler))

    # Wrong expected SHA-256 must raise ValueError, remove `.part`, and NOT overwrite with bad data
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        atomic_download(
            "https://example.org/weights.bin",
            dest,
            expected_sha256="0" * 64,
            client=client,
        )
    assert not dest.with_suffix(".bin.part").exists()

    # Correct SHA-256 detects existing `dest` is corrupt, re-downloads to `.part`, and atomically replaces!
    actual = atomic_download(
        "https://example.org/weights.bin",
        dest,
        expected_sha256=good_sha,
        expected_size=len(payload),
        client=client,
    )
    assert actual == good_sha
    assert dest.read_bytes() == payload
    assert not dest.with_suffix(".bin.part").exists()
