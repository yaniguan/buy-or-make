"""Tests for cost_aware_rank ordering, mixed-unit tiers, licenses, desalting, dedup, and filter_proposals."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import httpx

from buyormake import cost_aware_rank, evaluate_candidates, filter_proposals
from buyormake.canonical import normalize_smiles


def _write_stock_tsv(path: Path) -> None:
    rows = [
        "smiles\tvendor\tcatalog_id\tpack_value\tpack_unit\tprice\tcurrency\tsource_id\tsource_url\tlicense",
        # Buyable expensive: $50 / 1 g = $50/g
        "CC(=O)Oc1ccccc1C(=O)O\tVendorExp\tCAT-EXP\t1\tg\t50.0\tUSD\tlocal_stock\thttps://example.org/exp\tCC0-1.0",
        # Buyable cheap: $10 / 5 g = $2/g
        "CC(=O)Nc1ccc(O)cc1\tVendorCheap\tCAT-CHEAP\t5\tg\t10.0\tUSD\tlocal_stock\thttps://example.org/cheap\tCC0-1.0",
        # Buyable unknown price
        "c1ccccc1\tVendorNoPrice\tCAT-NP\t\t\t\t\tlocal_stock\thttps://example.org/np\tCC0-1.0",
        # Building blocks for 1-step amide: benzoic acid ($5/10g = $0.5/g) + aniline ($2/10g = $0.2/g)
        "O=C(O)c1ccccc1\tVendorBB\tBB-ACID\t10\tg\t5.0\tUSD\tlocal_stock\thttps://example.org/acid\tCC0-1.0",
        "Nc1ccccc1\tVendorBB\tBB-AMINE\t10\tg\t2.0\tUSD\tlocal_stock\thttps://example.org/amine\tCC0-1.0",
        # Building blocks for 2-step biaryl amide: 4-methoxybenzoic acid + 4-bromoaniline + pyridin-2-ylboronic acid (no price on boronic acid)
        "COc1ccc(C(=O)O)cc1\tVendorBB\tBB-MBA\t10\tg\t4.0\tUSD\tlocal_stock\thttps://example.org/mba\tCC0-1.0",
        "Nc1ccc(Br)cc1\tVendorBB\tBB-BRA\t10\tg\t8.0\tUSD\tlocal_stock\thttps://example.org/bra\tCC0-1.0",
        "OB(O)c1ccccn1\tVendorBB\tBB-PYB\t\t\t\t\tlocal_stock\thttps://example.org/pyb\tCC0-1.0",
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_ranking_order_and_invalid_smiles(tmp_path: Path) -> None:
    stock = tmp_path / "stock.tsv"
    _write_stock_tsv(stock)

    candidates = [
        "not_a_smiles",                           # invalid
        "C12C3C4C1C5C2C3C45",                     # cubane -> unknown (not in stock, no template route)
        "COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1",  # makeable in 2 steps (incomplete price)
        "O=C(Nc1ccccc1)c1ccccc1",                 # makeable in 1 step (complete price)
        "c1ccccc1",                               # buyable, price=None
        "CC(=O)Oc1ccccc1C(=O)O",                  # buyable, $50/g
        "CC(=O)Nc1ccc(O)cc1",                     # buyable, $2/g
    ]

    ranked = cost_aware_rank(
        candidates,
        max_steps=2,
        backends=("local_stock",),
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )

    statuses = [r.status for r in ranked]
    assert statuses == [
        "buyable",   # paracetamol $2/g
        "buyable",   # aspirin $50/g
        "buyable",   # benzene price=None
        "makeable",  # benzanilide 1-step
        "makeable",  # biaryl amide 2-step
        "unknown",   # cubane
        "invalid",   # not_a_smiles
    ]
    assert ranked[0].candidate == "CC(=O)Nc1ccc(O)cc1"
    # Priced buyable: price lives in best_price; there is no route cost estimate.
    assert ranked[0].best_price is not None
    assert ranked[0].est_cost is None
    assert ranked[0].cost_complete is False
    assert ranked[1].candidate == "CC(=O)Oc1ccccc1C(=O)O"
    assert ranked[2].candidate == "c1ccccc1"
    assert ranked[2].best_price is None
    assert ranked[2].cost_complete is False
    assert ranked[3].n_steps == 1
    assert ranked[3].cost_complete is True
    assert ranked[4].n_steps == 2
    assert ranked[4].cost_complete is False
    # The unpriced boronic acid leaf is named as excluded from the route cost.
    assert any(
        "leaves excluded from route cost" in r and "no listed USD/g price" in r
        for r in ranked[4].reasons
    )
    assert not any("leaves excluded from route cost" in r for r in ranked[3].reasons)
    assert ranked[6].inchikey is None


def test_rank_sort_key_never_mixes_usd_per_g_with_raw_pack_or_foreign_currency(tmp_path: Path) -> None:
    stock = tmp_path / "mixed_units.tsv"
    rows = [
        "smiles\tvendor\tcatalog_id\tpack_value\tpack_unit\tprice\tcurrency\tsource_id\tsource_url\tlicense",
        # Raw non-mass pack: $5 / 1 mL -> unit="USD" (value=5.0, NOT normalized USD/g!)
        "CCO\tV1\tC-RAW\t1\tmL\t5.0\tUSD\tlocal_stock\thttps://example.org/raw\tCC0-1.0",
        # Foreign currency: 10 EUR / 1 g -> unit="EUR/g" (value=10.0)
        "CCCO\tV1\tC-EUR\t1\tg\t10.0\tEUR\tlocal_stock\thttps://example.org/eur\tCC0-1.0",
        # Normalized USD/g: $12 / 1 g -> unit="USD/g" (value=12.0, MUST rank first!)
        "CCCCO\tV1\tC-USDG\t1\tg\t12.0\tUSD\tlocal_stock\thttps://example.org/usdg\tCC0-1.0",
        # Null price
        "c1ccccc1\tV1\tC-NULL\t\t\t\t\tlocal_stock\thttps://example.org/null\tCC0-1.0",
    ]
    stock.write_text("\n".join(rows) + "\n", encoding="utf-8")

    ranked = cost_aware_rank(
        ["CCO", "CCCO", "CCCCO", "c1ccccc1"],
        max_steps=1,
        backends=("local_stock",),
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    # Tier 0: USD/g ("CCCCO" at 12.0 USD/g), Tier 1: other units grouped by unit ("EUR/g", "USD"), Tier 2: null
    assert [r.candidate for r in ranked] == ["CCCCO", "CCCO", "CCO", "c1ccccc1"]
    assert ranked[0].best_price is not None and ranked[0].best_price.unit == "USD/g"
    assert ranked[1].best_price is not None and ranked[1].best_price.unit == "EUR/g"
    assert ranked[2].best_price is not None and ranked[2].best_price.unit == "USD"
    assert ranked[3].best_price is None


def test_filter_proposals_budget(tmp_path: Path) -> None:
    stock = tmp_path / "stock.tsv"
    _write_stock_tsv(stock)

    summary = filter_proposals(
        ["CC(=O)Nc1ccc(O)cc1", "CC(=O)Oc1ccccc1C(=O)O", "bad_smiles", "C12C3C4C1C5C2C3C45"],
        max_steps=2,
        budget=10.0,
        backends=("local_stock",),
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert [r.candidate for r in summary.kept] == ["CC(=O)Nc1ccc(O)cc1"]
    assert [r.candidate for r in summary.over_budget] == ["CC(=O)Oc1ccccc1C(=O)O"]
    assert len(summary.rejected) == 3
    over_dec = next(d for d in summary.decisions if d.candidate == "CC(=O)Oc1ccccc1C(=O)O")
    assert over_dec.over_budget is True
    assert over_dec.keep is False
    assert any("over_budget" in r for r in over_dec.reasons)


def test_filter_proposals_unit_currency_and_completeness_policies(tmp_path: Path) -> None:
    stock = tmp_path / "policy_stock.tsv"
    rows = [
        "smiles\tvendor\tcatalog_id\tpack_value\tpack_unit\tprice\tcurrency\tsource_id\tsource_url\tlicense",
        # 500 JPY/g (~$3.30/g, NOT > 25 USD/g, but incompatible with USD budget!)
        "CCO\tVJPY\tJ-1\t1\tg\t500.0\tJPY\tlocal_stock\thttps://example.org/j1\tCC0-1.0",
        # $20 / 1 mL vial -> unit="USD" raw pack (NOT comparable to per-gram USD/g budget!)
        "CCCO\tVRaw\tR-1\t1\tmL\t20.0\tUSD\tlocal_stock\thttps://example.org/r1\tCC0-1.0",
        # Unpriced buyable
        "c1ccccc1\tVNull\tN-1\t\t\t\t\tlocal_stock\thttps://example.org/n1\tCC0-1.0",
        # Valid $8/g USD offer
        "CCCCO\tVOK\tOK-1\t1\tg\t8.0\tUSD\tlocal_stock\thttps://example.org/ok1\tCC0-1.0",
    ]
    stock.write_text("\n".join(rows) + "\n", encoding="utf-8")

    # Strict budget gate (`require_complete_cost=True`, `budget_currency="USD"`)
    strict = filter_proposals(
        ["CCO", "CCCO", "c1ccccc1", "CCCCO"],
        budget=25.0,
        budget_currency="USD",
        require_complete_cost=True,
        backends=("local_stock",),
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert [r.candidate for r in strict.kept] == ["CCCCO"]
    # Neither 500 JPY/g nor $20 raw-pack nor unpriced is falsely marked `over_budget`
    assert strict.over_budget == []
    by_cand = {d.candidate: d for d in strict.decisions}
    assert by_cand["CCO"].price_status == "incompatible_unit"
    assert by_cand["CCO"].over_budget is False
    assert by_cand["CCO"].keep is False
    assert by_cand["CCCO"].price_status == "incompatible_unit"
    assert by_cand["CCCO"].keep is False
    assert by_cand["c1ccccc1"].price_status == "missing"
    assert by_cand["c1ccccc1"].keep is False

    # Same 500 JPY/g offer passes under JPY budget=1000.0
    jpy_summary = filter_proposals(
        ["CCO"],
        budget=1000.0,
        budget_currency="JPY",
        require_complete_cost=True,
        backends=("local_stock",),
        offline=True,
        cache_path=tmp_path / "cache_jpy.sqlite",
        stock_paths=[stock],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert [r.candidate for r in jpy_summary.kept] == ["CCO"]


def test_route_and_leaf_licenses_propagate_into_availability_record(tmp_path: Path) -> None:
    target = "CC(=O)Nc1ccc(F)cc1"
    leaf1 = "CC(=O)O"
    leaf2 = "Nc1ccc(F)cc1"
    t_nm = normalize_smiles(target)
    l1_nm = normalize_smiles(leaf1)
    l2_nm = normalize_smiles(leaf2)
    assert t_nm and l1_nm and l2_nm

    rxn_db = tmp_path / "ord_rxn.sqlite"
    conn = sqlite3.connect(str(rxn_db))
    conn.execute(
        "CREATE TABLE reactions (id INTEGER PRIMARY KEY, product_inchikey TEXT, "
        "product_smiles TEXT, reactants_json TEXT, reaction_smiles TEXT, reference TEXT, source TEXT)"
    )
    conn.execute(
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, reaction_smiles, reference, source) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        (
            t_nm.inchikey,
            t_nm.canonical_smiles,
            json.dumps([leaf1, leaf2]),
            f"{leaf1}.{leaf2}>>{target}",
            "ord-rxn-001",
            "ord_lookup",
        ),
    )
    conn.commit()
    conn.close()

    def handler(request: httpx.Request) -> httpx.Response:
        url = str(request.url)
        # Target not in PubChem; leaves are in PubChem
        if l1_nm.inchikey in url:
            return httpx.Response(200, json={"IdentifierList": {"CID": [176]}})
        if l2_nm.inchikey in url:
            return httpx.Response(200, json={"IdentifierList": {"CID": [9794]}})
        if "/categories/compound/176/JSON" in url or "/categories/compound/9794/JSON" in url:
            return httpx.Response(
                200,
                json={
                    "SourceCategories": {
                        "Categories": [
                            {
                                "Category": "Chemical Vendors",
                                "Sources": [{"SourceName": "LeafVendor", "RegistryID": "LV-1"}],
                            }
                        ]
                    }
                },
            )
        return httpx.Response(404)

    ranked, records = evaluate_candidates(
        [target],
        max_steps=1,
        backends=("pubchem",),
        route_backend="lookup",
        offline=False,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[],
        reaction_db_path=rxn_db,
        scscore_weights=tmp_path / "no_sc.json.gz",
        http_client=httpx.Client(transport=httpx.MockTransport(handler)),
    )
    assert ranked[0].status == "makeable"
    rec = records[0]
    assert "ord_data" in rec.provenance.used
    assert "pubchem_pug" in rec.provenance.used
    assert "CC-BY-SA-4.0" in rec.licenses
    assert "custom/terms" in rec.licenses


def test_desalting_and_candidate_deduplication(tmp_path: Path, monkeypatch) -> None:
    stock = tmp_path / "lidocaine.smi"
    # Neutral lidocaine only in stock
    stock.write_text("CCN(CC)CC(=O)Nc1c(C)cccc1C lidocaine\nCCO ethanol\n", encoding="utf-8")

    # Lidocaine HCl salt matches neutral lidocaine parent; duplicate ethanol SMILES reuse evaluation
    lidocaine_hcl = "CCN(CC)CC(=O)Nc1c(C)cccc1C.Cl"
    parent_nm = normalize_smiles(lidocaine_hcl, desalt=True)
    assert parent_nm is not None
    assert parent_nm.canonical_smiles == "CCN(CC)CC(=O)Nc1c(C)cccc1C"

    import buyormake.rank as rank_mod

    score_calls = {"n": 0}
    orig_compute_scores = rank_mod.compute_scores

    def counted_scores(*args, **kwargs):
        score_calls["n"] += 1
        return orig_compute_scores(*args, **kwargs)

    monkeypatch.setattr(rank_mod, "compute_scores", counted_scores)

    ranked, records = evaluate_candidates(
        [lidocaine_hcl, "CCO", "OCC", "C(C)O"],
        candidate_names=["LidoHCl", "Et1", "Et2", "Et3"],
        max_steps=1,
        backends=("local_stock",),
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert len(ranked) == 4
    assert len(records) == 4
    # Only 2 unique InChIKeys (lidocaine HCl and ethanol) -> compute_scores called twice, not 4 times!
    assert score_calls["n"] == 2
    lido_res = next(r for r in ranked if r.candidate == lidocaine_hcl)
    assert lido_res.status == "buyable"
    assert any("desalted parent structure" in r for r in lido_res.reasons)
