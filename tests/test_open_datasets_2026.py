"""Tests for 2024-2026 open chemistry dataset registry, SPDX invariants, and multi-format parsers."""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sqlite3
import zipfile
from pathlib import Path

import httpx
import pytest
from pydantic import ValidationError

from buyormake import (
    CORPUS_ARTIFACTS,
    OPEN_DATASET_2026_IDS,
    SOURCE_BY_ID,
    SOURCES,
    Source,
    download_corpus_artifact,
    evaluate_candidates,
    import_chembl_chemreps,
    import_coconut_compounds,
    import_paroutes_stock,
    import_surechembl_compounds,
    import_synplanner_building_blocks,
    import_syntheseus_inventory,
    index_orderly_reactions,
    index_paroutes_reactions_csv,
    index_paroutes_routes,
    index_synplanner_reactions,
    index_uspto_50k_reactions,
    load_paroutes_routes,
    load_synplanner_targets,
    parse_chembl_chemreps_records,
    parse_coconut_records,
    parse_paroutes_route_tree,
    select_record_license,
)
from buyormake.canonical import normalize_smiles
from buyormake.open_corpora import ingest_downloaded_artifact
from buyormake.routes import ReactionLookupBackend
from buyormake.stock import LocalStockBackend

REPO_ROOT = Path(__file__).resolve().parents[1]


def _has_cjk(text: str) -> bool:
    return any(
        "\u4e00" <= ch <= "\u9fff"
        or "\u3400" <= ch <= "\u4dbf"
        or "\u3040" <= ch <= "\u30ff"
        or "\uac00" <= ch <= "\ud7af"
        for ch in text
    )


def test_open_datasets_2026_source_registry_and_spdx_invariants() -> None:
    expected_spdx = {
        "paroutes": "CC-BY-4.0",
        "chembl": "CC-BY-SA-3.0",
        "surechembl": "CC-BY-4.0",
        "coconut": "CC0-1.0",
        "orderly": "CC-BY-4.0",
        "syntheseus_uspto_50k": "MIT",
        "synplanner_data": "MIT",
    }
    assert set(OPEN_DATASET_2026_IDS) == set(expected_spdx.keys())
    assert len(OPEN_DATASET_2026_IDS) >= 5
    assert len(SOURCES) >= 19

    for sid, spdx in expected_spdx.items():
        assert sid in SOURCE_BY_ID
        src = SOURCE_BY_ID[sid]
        assert src.spdx_license == spdx
        assert src.license == spdx
        assert src.license_checked == "2026-11-25"
        assert src.license_check_method == "direct"
        assert src.commercial_use_ok is True
        assert src.redistribution_ok is True
        assert src.status == "fetch_only"
        assert src.homepage.startswith("https://")
        assert src.license_evidence.startswith("https://")
        assert src.license_quote is not None and len(src.license_quote) >= 40

    # All 19 sources carry populated non-empty spdx_license and zero CJK characters
    for src in SOURCES:
        assert src.spdx_license is not None and len(src.spdx_license.strip()) > 0
        dumped = src.model_dump_json()
        assert not _has_cjk(dumped), f"CJK character found in Source {src.id}"

    # Source.spdx_license default population and blank rejection
    auto_mit = Source(
        id="auto_mit",
        name="Auto MIT",
        homepage="https://example.org",
        license="MIT",
        license_evidence="https://example.org/license",
        license_checked="2026-11-25",
        commercial_use_ok=True,
        redistribution_ok=True,
        attribution="Example",
        access="github",
    )
    assert auto_mit.spdx_license == "MIT"

    auto_custom = Source(
        id="auto_custom",
        name="Auto Custom",
        homepage="https://example.org",
        license="custom/terms",
        license_evidence="https://example.org/terms",
        license_checked="2026-11-25",
        commercial_use_ok=False,
        redistribution_ok=False,
        attribution="Example",
        access="download",
    )
    assert auto_custom.spdx_license == "LicenseRef-Custom-Terms"

    with pytest.raises(ValidationError):
        Source(
            id="blank_spdx",
            name="Blank SPDX",
            homepage="https://example.org",
            license="MIT",
            spdx_license="   ",
            license_evidence="https://example.org/license",
            license_checked="2026-11-25",
            commercial_use_ok=True,
            redistribution_ok=True,
            attribution="Example",
            access="github",
        )

    # Governing license restrictiveness handles CC-BY-SA-3.0 and LicenseRef-Custom-Terms
    assert select_record_license(["CC0-1.0", "MIT", "CC-BY-4.0"]) == "CC-BY-4.0"
    assert select_record_license(["CC-BY-4.0", "CC-BY-SA-3.0"]) == "CC-BY-SA-3.0"
    assert select_record_license(["CC-BY-SA-3.0", "LicenseRef-Custom-Terms"]) == "LicenseRef-Custom-Terms"


def test_corpus_artifacts_registry_and_atomic_download(tmp_path: Path, monkeypatch) -> None:
    covered_sources = {art.source_id for art in CORPUS_ARTIFACTS.values()}
    assert set(OPEN_DATASET_2026_IDS).issubset(covered_sources)
    for art in CORPUS_ARTIFACTS.values():
        assert art.url.startswith("https://")
        assert art.source_id in SOURCE_BY_ID

    payload = b"CCO\nc1ccccc1\n"
    md5_hex = hashlib.md5(payload).hexdigest()
    sha_hex = hashlib.sha256(payload).hexdigest()

    from buyormake.open_corpora import CorpusArtifact

    test_spec = CorpusArtifact(
        artifact_key="test_paroutes_mini",
        source_id="paroutes",
        description="Mini test stock",
        url="https://zenodo.org/api/records/7341155/files/stock_n1.txt/content",
        filename="stock_n1.txt",
        kind="stock",
        expected_md5=md5_hex,
        expected_size=len(payload),
    )
    monkeypatch.setitem(CORPUS_ARTIFACTS, "test_paroutes_mini", test_spec)

    client = httpx.Client(
        transport=httpx.MockTransport(
            lambda req: httpx.Response(200, content=payload, headers={"Content-Length": str(len(payload))})
        )
    )
    out_dir = tmp_path / "data" / "open_corpora"
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        download_corpus_artifact(
            "test_paroutes_mini",
            out_root=out_dir,
            expected_sha256="0" * 64,
            client=client,
        )
    assert not (out_dir / "paroutes" / "stock_n1.txt").exists()

    dest = download_corpus_artifact(
        "test_paroutes_mini",
        out_root=out_dir,
        expected_sha256=sha_hex,
        client=client,
    )
    assert dest == out_dir / "paroutes" / "stock_n1.txt"
    assert dest.read_bytes() == payload

    with pytest.raises(KeyError):
        download_corpus_artifact("nonexistent_artifact_key", out_root=out_dir, client=client)


def test_paroutes_stock_and_multistep_route_tree_triage(tmp_path: Path) -> None:
    # Target: N-(4-(pyridin-2-yl)phenyl)acetamide via 2-step PaRoutes tree:
    # Step 1 (amide): CC(=O)Nc1ccc(-c2ccccn2)cc1 <- CC(=O)O + Nc1ccc(-c2ccccn2)cc1
    # Step 2 (Suzuki): Nc1ccc(-c2ccccn2)cc1 <- Nc1ccc(Br)cc1.Cl (salt in stock!) + OB(O)c1ccccn1
    stock_txt = tmp_path / "stock_n1.txt"
    stock_txt.write_text(
        "# PaRoutes n1 stock fixture\n"
        "CC(=O)O\n"
        "Nc1ccc(Br)cc1.Cl\n"
        "OB(O)c1ccccn1\n",
        encoding="utf-8",
    )
    stock_db = tmp_path / "stock.sqlite"
    n_stock = import_paroutes_stock(stock_txt, db_path=stock_db, subset="n1")
    assert n_stock == 3

    target_smi = "CC(=O)Nc1ccc(-c2ccccn2)cc1"
    inter_smi = "Nc1ccc(-c2ccccn2)cc1"
    route_tree = {
        "type": "mol",
        "smiles": target_smi,
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "metadata": {
                    "ID": "US20240001A1",
                    "rsmi": f"CC(=O)O.{inter_smi}.CN(C)C=O>>{target_smi}",
                },
                "children": [
                    {"type": "mol", "smiles": "CC(=O)O", "in_stock": True},
                    {
                        "type": "mol",
                        "smiles": inter_smi,
                        "in_stock": False,
                        "children": [
                            {
                                "type": "reaction",
                                "metadata": {
                                    "ID": "US20240002A1",
                                    "rsmi": f"Nc1ccc(Br)cc1.OB(O)c1ccccn1.[Pd]>>{inter_smi}",
                                },
                                "children": [
                                    {"type": "mol", "smiles": "Nc1ccc(Br)cc1", "in_stock": True},
                                    {"type": "mol", "smiles": "OB(O)c1ccccn1", "in_stock": True},
                                ],
                            }
                        ],
                    },
                ],
            }
        ],
    }

    parsed_route = parse_paroutes_route_tree(route_tree)
    assert parsed_route is not None
    assert parsed_route.n_steps == 2
    assert parsed_route.leaves_purchasable == [True, True, True]

    routes_gz = tmp_path / "ref_routes_n1.json.gz"
    with gzip.open(routes_gz, "wt", encoding="utf-8") as gz:
        json.dump([route_tree], gz)

    loaded_routes = load_paroutes_routes(routes_gz)
    assert len(loaded_routes) == 1

    rxn_db = tmp_path / "reactions.sqlite"
    n_indexed = index_paroutes_routes(routes_gz, db_path=rxn_db)
    assert n_indexed == 2

    # Also test selected_reactions.csv indexing
    rxn_csv = tmp_path / "selected_reactions.csv"
    rxn_csv.write_text(
        "id,rxn_smiles\n"
        "paroutes_csv_1,CC(=O)Cl.Nc1ccccc1.CCN(CC)CC>>CC(=O)Nc1ccccc1\n",
        encoding="utf-8",
    )
    assert index_paroutes_reactions_csv(rxn_csv, db_path=rxn_db) == 1

    ranked, records = evaluate_candidates(
        [target_smi],
        max_steps=3,
        backends=("local_stock",),
        route_backend="lookup",
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock_db],
        reaction_db_path=rxn_db,
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert len(ranked) == 1
    assert ranked[0].status == "makeable"
    assert ranked[0].n_steps == 2
    assert ranked[0].est_cost is None
    assert "paroutes" in records[0].provenance.used
    assert "CC-BY-4.0" in records[0].licenses
    af_rec = records[0].to_assayfold_record()
    assert af_rec.license == "CC-BY-4.0"
    assert af_rec.record_hash == af_rec.content_hash()


def test_chembl_surechembl_and_coconut_importers_and_assayfold_sealing(
    tmp_path: Path, monkeypatch
) -> None:
    # 1. ChEMBL 37 chemreps .txt.gz
    chembl_gz = tmp_path / "chembl_37_chemreps.txt.gz"
    with gzip.open(chembl_gz, "wt", encoding="utf-8") as gz:
        gz.write(
            "chembl_id\tcanonical_smiles\tstandard_inchi\tstandard_inchi_key\n"
            "CHEMBL25\tCC(=O)Oc1ccccc1C(=O)O\tInChI=1S/...\tBSYNRYMUTXBXSQ-UHFFFAOYSA-N\n"
            "CHEMBL545\tCCO\tInChI=1S/...\tLFQSCWFLJHTTHZ-UHFFFAOYSA-N\n"
        )

    ref_db = tmp_path / "reference_compounds.sqlite"
    stock_db = tmp_path / "stock.sqlite"
    assert import_chembl_chemreps(chembl_gz, db_path=ref_db) == 2
    assert not stock_db.exists(), "Reference compound importers must not touch stock.sqlite by default"

    chembl_records = parse_chembl_chemreps_records(chembl_gz)
    assert len(chembl_records) == 2
    assert chembl_records[0].license == "CC-BY-SA-3.0"
    assert chembl_records[0].source_id == "chembl"
    assert chembl_records[0].record_hash == chembl_records[0].content_hash()

    # 2. SureChEMBL TSV + mocked Parquet path
    sc_tsv = tmp_path / "surechembl_sample.tsv"
    sc_tsv.write_text(
        "SureChEMBL_ID\tSMILES\n"
        "SCHEMBL1001\tc1ccc(-c2ccccn2)cc1\n",
        encoding="utf-8",
    )
    assert import_surechembl_compounds(sc_tsv, db_path=ref_db) == 1

    fake_pq = tmp_path / "compounds.parquet"
    fake_pq.write_bytes(b"parquet_stub")
    monkeypatch.setattr(
        "buyormake.open_corpora._read_parquet_columns",
        lambda _: {"surechembl_id": ["SCHEMBL2002"], "smiles": ["COc1ccccc1"]},
    )
    assert import_surechembl_compounds(fake_pq, db_path=ref_db) == 1

    # 3. COCONUT 2.0 .zip archive containing CSV
    coco_zip = tmp_path / "coconut_csv_lite-11-2026.zip"
    with zipfile.ZipFile(coco_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            "coconut_csv_lite.csv",
            "identifier,canonical_smiles,name\n"
            "CNP0100001.1,COc1cc(/C=C/C(=O)O)ccc1O,Ferulic acid\n"
            "CNP0100002.1,Cn1cnc2c1c(=O)n(C)c(=O)n2C,Caffeine\n",
        )
    assert import_coconut_compounds(coco_zip, db_path=ref_db) == 2

    coco_records = parse_coconut_records(coco_zip)
    assert len(coco_records) == 2
    assert coco_records[0].license == "CC0-1.0"
    assert coco_records[0].payload.name == "Ferulic acid"
    assert coco_records[0].record_hash == coco_records[0].content_hash()

    # Verify LocalStockBackend reads reference_compounds.sqlite with unpriced honesty & SPDX licenses
    backend = LocalStockBackend([ref_db])
    try:
        caff_nm = normalize_smiles("Cn1cnc2c1c(=O)n(C)c(=O)n2C")
        asp_nm = normalize_smiles("CC(=O)Oc1ccccc1C(=O)O")
        sc_nm = normalize_smiles("COc1ccccc1")
        assert caff_nm and asp_nm and sc_nm
        caff_offs = backend.lookup(caff_nm.inchikey)
        asp_offs = backend.lookup(asp_nm.inchikey)
        sc_offs = backend.lookup(sc_nm.inchikey)
        assert len(caff_offs) == 1 and caff_offs[0].license == "CC0-1.0" and caff_offs[0].price is None
        assert len(asp_offs) == 1 and asp_offs[0].license == "CC-BY-SA-3.0" and asp_offs[0].in_stock is None
        assert len(sc_offs) == 1 and sc_offs[0].license == "CC-BY-4.0" and sc_offs[0].source_id == "surechembl"
    finally:
        backend.close()


def test_orderly_syntheseus_and_synplanner_reaction_and_stock_adapters(tmp_path: Path) -> None:
    rxn_db = tmp_path / "reactions.sqlite"
    stock_db = tmp_path / "stock.sqlite"

    # 1. ORDerly CSV with separate reactants, products, and unmapped DMF solvent spectator
    orderly_csv = tmp_path / "orderly_retro_test.csv"
    orderly_csv.write_text(
        "reaction_id,reactants,products,solvents\n"
        "ord_clean_01,[CH3:1][C:2](=[O:3])[OH:9].[NH2:4][c:5]1[cH:6][cH:7][c:8]([Cl:10])[cH:11][cH:12]1,"
        "[CH3:1][C:2](=[O:3])[NH:4][c:5]1[cH:6][cH:7][c:8]([Cl:10])[cH:11][cH:12]1,CN(C)C=O\n",
        encoding="utf-8",
    )
    assert index_orderly_reactions(orderly_csv, db_path=rxn_db) == 1

    # 2. Syntheseus USPTO-50K .zip (3-part reactants>reagents>products) + Retro* .csv.gz inventory
    uspto50k_zip = tmp_path / "uspto50k.zip"
    with zipfile.ZipFile(uspto50k_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            "test.csv",
            "id,class,reactions\n"
            "US50K_001,1,OB(O)c1ccccc1.Brc1ccc(C#N)cc1>C1CCOC1.[Pd]>N#Cc1ccc(-c2ccccc2)cc1\n",
        )
    assert index_uspto_50k_reactions(uspto50k_zip, db_path=rxn_db) == 1

    inv_gz = tmp_path / "ORIGINAL-inventory.csv.gz"
    with gzip.open(inv_gz, "wt", encoding="utf-8") as gz:
        gz.write(
            "id,smiles\n"
            "INV-1,OB(O)c1ccccc1\n"
            "INV-2,Brc1ccc(C#N)cc1\n"
        )
    assert import_syntheseus_inventory(inv_gz, db_path=stock_db) == 2

    # 3. SynPlanner building blocks (filtering bare single-atom ions), targets, and standardized USPTO zip
    sp_bb = tmp_path / "building_blocks_em_sa_ln.smi"
    sp_bb.write_text(
        "[Li+]\n"
        "[Na+]\n"
        "Cl\n"
        "CC(=O)O\n"
        "Nc1ccc(Cl)cc1\n",
        encoding="utf-8",
    )
    assert import_synplanner_building_blocks(sp_bb, db_path=stock_db) == 2

    sp_targets = tmp_path / "targets_with_sascore_1.5_2.5.smi"
    sp_targets.write_text("N#Cc1ccc(-c2ccccc2)cc1\nCC(=O)Nc1ccc(Cl)cc1\n", encoding="utf-8")
    targets = load_synplanner_targets(sp_targets)
    assert len(targets) == 2

    sp_rxn_zip = tmp_path / "uspto_standardized.zip"
    with zipfile.ZipFile(sp_rxn_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(
            "uspto_std.smi",
            "CC(=O)O.Nc1ccccc1>>CC(=O)Nc1ccccc1 US_SP_1\n",
        )
    assert index_synplanner_reactions(sp_rxn_zip, db_path=rxn_db) == 1

    # Verify ReactionLookupBackend stripped spectators and preserved stoichiometric reactants
    lookup = ReactionLookupBackend(rxn_db)
    try:
        p_biaryl = normalize_smiles("N#Cc1ccc(-c2ccccc2)cc1")
        assert p_biaryl is not None
        steps = lookup.one_step(p_biaryl.inchikey, p_biaryl.canonical_smiles)
        assert len(steps) == 1
        assert set(steps[0].reactants) == {"OB(O)c1ccccc1", "N#Cc1ccc(Br)cc1"}
        assert steps[0].source == "syntheseus_uspto_50k_lookup"
    finally:
        lookup.close()

    conn = sqlite3.connect(str(rxn_db))
    agents = [r[0] for r in conn.execute("SELECT agents_json FROM reactions").fetchall()]
    conn.close()
    assert any("CN(C)C=O" in a for a in agents)
    assert any("C1CCOC1" in a for a in agents)

    # End-to-end offline triage over SynPlanner targets using Syntheseus + ORDerly reactions & stock
    ranked, records = evaluate_candidates(
        targets,
        max_steps=2,
        backends=("local_stock",),
        route_backend="lookup",
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock_db],
        reaction_db_path=rxn_db,
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert [r.status for r in ranked] == ["makeable", "makeable"]
    used_union = set().union(*(r.provenance.used for r in records))
    assert {"orderly", "syntheseus_uspto_50k", "synplanner_data"}.issubset(used_union)

    # Verify ingest_downloaded_artifact & CLI fetcher script smoke path
    assert ingest_downloaded_artifact("synplanner_targets_sa_1_5_2_5", sp_targets) == 2

    script_path = REPO_ROOT / "scripts" / "fetch_open_datasets_2026.py"
    spec = importlib.util.spec_from_file_location("fetch_open_datasets_2026", script_path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.main(["--list"]) == 0
    cli_rxn_db = tmp_path / "cli_rxn.sqlite"
    assert (
        mod.main(
            [
                "--artifact",
                "syntheseus_uspto_50k",
                "--local-file",
                str(uspto50k_zip),
                "--reaction-db",
                str(cli_rxn_db),
            ]
        )
        == 0
    )
    assert cli_rxn_db.exists()
