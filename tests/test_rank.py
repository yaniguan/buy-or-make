"""Tests for cost_aware_rank ordering, invalid SMILES, and filter_proposals."""

from __future__ import annotations

from pathlib import Path

from buyormake import cost_aware_rank, filter_proposals


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
