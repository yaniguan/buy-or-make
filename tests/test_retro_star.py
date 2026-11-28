"""Tests for `Retro*` best-first AND-OR tree search, cycle prevention, extended SMARTS, and discrete MOQ."""

from __future__ import annotations

import json
import sqlite3
from collections import Counter
from pathlib import Path

import pytest
from rdkit import Chem

from buyormake import (
    Offer,
    Quantity,
    RetroStarPlanner,
    evaluate_candidates,
    optimize_discrete_pack_order,
    optimize_route_moq_outlay,
)
from buyormake.canonical import normalize_smiles
from buyormake.cli import main as cli_main
from buyormake.retro_star import (
    EXTENDED_RETRO_STAR_TEMPLATES,
    RETRO_STAR_TEMPLATE_DELTAS,
)
from buyormake.schema import RouteStep


def _heavy_counts(smi: str) -> Counter[str]:
    mol = Chem.MolFromSmiles(smi)
    assert mol is not None, f"Unparseable SMILES: {smi}"
    return Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def test_all_extended_retro_star_templates_conserve_exact_element_deltas() -> None:
    from rdkit.Chem import rdChemReactions

    planner = RetroStarPlanner(max_branch_per_node=60)
    assert len(EXTENDED_RETRO_STAR_TEMPLATES) == 100
    seen_names: set[str] = set()
    for tpl_name, smarts in EXTENDED_RETRO_STAR_TEMPLATES:
        assert tpl_name not in seen_names, f"Duplicate template name: {tpl_name}"
        seen_names.add(tpl_name)
        assert tpl_name in RETRO_STAR_TEMPLATE_DELTAS, f"Missing delta for {tpl_name}"
        rxn = rdChemReactions.ReactionFromSmarts(smarts)
        assert rxn is not None, f"Failed to compile SMARTS for {tpl_name}"
        assert rxn.GetNumReactantTemplates() == 1

    probes = [
        "Fc1ccc(-c2ccccn2)cc1",
        "CC(=O)Nc1ccc(-c2ccc(C(=O)OC)cc2)cc1",
        "CC(C)(C)OC(=O)N1CCN(c2ncnc3ccccc23)CC1",
        "Nc1ccc(S(=O)(=O)N2CCOCC2)cc1",
        "O=C(O)c1ccc(OCc2ccccc2)cc1",
        "CC(O)c1ccc(C#Cc2ccccc2)cc1",
        "CCN(CCc1ccccc1)Cc1ccc(Cl)cc1",
        "CSCc1ccc(O)cc1",
        "O=C1CCCO1",
        "O=C(Nc1ccc(F)cc1)N1CCOCC1",
        "CCOC(=O)N1CCN(Cc2ccccc2)CC1",
        "CCOC(=O)/C=C/c1ccc(OC)cc1",
        "Cc1csc(Nc2ccccc2)n1",
        "Cc1nc2ccccc2[nH]1",
        "Cc1noc(-c2ccccc2)n1",
        "c1ccc(-c2cn(Cc3ccccc3)nn2)cc1",
        "OB(O)c1ccc(F)cc1",
        "O=C(Cl)c1ccc(Cl)cc1",
        "O=C=Nc1ccc(OC)cc1",
        "CC1CCN(Cc2ccccc2)CC1",
    ]
    total_checked = 0
    for probe in probes:
        p_counts = _heavy_counts(probe)
        for cost_val, step in planner.expand_extended_smarts(probe, allow_unstocked_fgi=True):
            assert cost_val > 0.0
            tpl_name = step.template_or_reference.split(":", 1)[0].strip()
            assert tpl_name in RETRO_STAR_TEMPLATE_DELTAS
            expected = Counter(p_counts)
            for sym, delta in RETRO_STAR_TEMPLATE_DELTAS[tpl_name].items():
                expected[sym] += delta
                if expected[sym] <= 0:
                    del expected[sym]
            actual = sum((_heavy_counts(r) for r in step.reactants), Counter())
            assert actual == expected, f"Delta mismatch on {tpl_name}"
            total_checked += 1
    assert total_checked >= 35


def test_retro_star_solves_convergent_and_deep_4step_routes_missed_by_greedy_2step(
    tmp_path: Path,
) -> None:
    # Build a 4-step CONVERGENT synthesis target where Step 1 splits T into TWO non-purchasable
    # intermediates (I_left and I_right), I_left requires 2 steps, and I_right requires 1 step:
    #   Step 1: T -> I_left + I_right (convergent!)
    #   Step 2: I_left -> M_sub + B1
    #   Step 3: M_sub -> B2 + B3
    #   Step 4: I_right -> B4 + B5
    target = "COc1ccc(C(=O)Nc2ccc(-c3ccc(C(=O)NCCc4ccccn4)cc3)cc2)cc1"
    i_left = "COc1ccc(C(=O)Nc2ccc(Br)cc2)cc1"
    i_right = "OB(O)c1ccc(C(=O)NCCc2ccccn2)cc1"
    m_sub = "COc1ccc(C(=O)Cl)cc1"
    b1 = "Nc1ccc(Br)cc1"
    b2 = "COc1ccc(C(=O)O)cc1"
    b3 = "O=S(Cl)Cl"
    b4 = "OB(O)c1ccc(C(=O)O)cc1"
    b5 = "NCCc1ccccn1"

    rxn_db = tmp_path / "convergent_4step.sqlite"
    conn = sqlite3.connect(str(rxn_db))
    conn.execute(
        "CREATE TABLE reactions (id INTEGER PRIMARY KEY, product_inchikey TEXT, "
        "product_smiles TEXT, reactants_json TEXT, reaction_smiles TEXT, reference TEXT, source TEXT)"
    )
    steps_spec = [
        (target, [i_left, i_right], "paroutes:conv1"),
        (i_left, [m_sub, b1], "paroutes:conv2"),
        (m_sub, [b2, b3], "paroutes:conv3"),
        (i_right, [b4, b5], "paroutes:conv4"),
    ]
    for prod_s, reacts, ref in steps_spec:
        p_nm = normalize_smiles(prod_s)
        assert p_nm is not None
        conn.execute(
            "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
            "reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?)",
            (
                p_nm.inchikey,
                p_nm.canonical_smiles,
                json.dumps(reacts),
                f"{'.'.join(reacts)}>>{p_nm.canonical_smiles}",
                ref,
                "paroutes_lookup",
            ),
        )
    # Also inject a directed cycle trap: M_sub <-> cycle_node to verify cycle prevention!
    m_nm = normalize_smiles(m_sub)
    cyc_smi = "COc1ccc(C(=O)Br)cc1"
    cyc_nm = normalize_smiles(cyc_smi)
    assert m_nm and cyc_nm
    conn.execute(
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
        "reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?)",
        (m_nm.inchikey, m_nm.canonical_smiles, json.dumps([cyc_smi]), f"{cyc_smi}>>{m_sub}", "cyc1", "paroutes_lookup"),
    )
    conn.execute(
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
        "reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?)",
        (cyc_nm.inchikey, cyc_nm.canonical_smiles, json.dumps([m_sub]), f"{m_sub}>>{cyc_smi}", "cyc2", "paroutes_lookup"),
    )
    conn.commit()
    conn.close()

    stock_tsv = tmp_path / "stock.tsv"
    stock_rows = [
        "smiles\tvendor\tcatalog_id\tpack_value\tpack_unit\tprice\tcurrency\tsource_id\tsource_url\tlicense",
    ]
    for idx, b_smi in enumerate([b1, b2, b3, b4, b5]):
        stock_rows.append(
            f"{b_smi}\tBBVendor\tBB-{idx}\t1.0\tg\t2.50\tUSD\tlocal_stock\thttps://ex/{idx}\tCC0-1.0"
        )
    stock_tsv.write_text("\n".join(stock_rows) + "\n", encoding="utf-8")

    # 1. Legacy lookup_smarts fails on convergent 4-step route (both Step-1 reactants out of stock, depth > 3)
    ranked_legacy, _ = evaluate_candidates(
        [target],
        max_steps=5,
        backends=("local_stock",),
        route_backend="lookup",
        planner="lookup_smarts",
        offline=True,
        cache_path=tmp_path / "c_leg.sqlite",
        stock_paths=[stock_tsv],
        reaction_db_path=rxn_db,
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert ranked_legacy[0].status == "unknown"

    # 2. Retro* solves the convergent 4-step AND-OR tree cleanly without getting trapped in the cycle!
    ranked_rs, recs_rs = evaluate_candidates(
        [target],
        max_steps=5,
        backends=("local_stock",),
        route_backend="lookup",
        planner="retro_star",
        offline=True,
        cache_path=tmp_path / "c_rs.sqlite",
        stock_paths=[stock_tsv],
        reaction_db_path=rxn_db,
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert ranked_rs[0].status == "makeable"
    assert ranked_rs[0].n_steps == 4
    assert ranked_rs[0].cost_complete is True
    assert ranked_rs[0].data_source == "local_vendor_catalog"
    best_route = recs_rs[0].routes[0]
    assert best_route.n_steps == 4
    assert best_route.tree_depth == 3  # Convergent: 4 total reactions, longest linear chain = 3!
    assert best_route.score is not None and 0.0 < best_route.score <= 1.0

    # 3. CLI `--planner retro_star --max-steps 5` end-to-end check
    smi_file = tmp_path / "query.smi"
    smi_file.write_text(f"{target} ConvTarget\n", encoding="utf-8")
    out_jsonl = tmp_path / "out.jsonl"
    rc = cli_main(
        [
            "rank",
            str(smi_file),
            "--planner",
            "retro_star",
            "--route-backend",
            "lookup",
            "--max-steps",
            "5",
            "--backends",
            "local_stock",
            "--offline",
            "--stock",
            str(stock_tsv),
            "--reaction-db",
            str(rxn_db),
            "--records-out",
            str(out_jsonl),
            "--json",
        ]
    )
    assert rc == 0
    assert out_jsonl.exists()


def test_exact_mixed_sku_and_shared_synthon_discrete_moq_optimization() -> None:
    acid = "O=C(O)c1ccccc1"
    acid_nm = normalize_smiles(acid)
    assert acid_nm is not None

    offers = [
        Offer(
            inchikey=acid_nm.inchikey,
            vendor="ChemCo",
            catalog_id="VIAL-250MG",
            pack_size=Quantity(value=250.0, unit="mg", raw_unit_str="250 mg"),
            price=18.0,
            currency="USD",
            price_basis="listed_catalog",
            source_id="local_stock",
            source_url="https://ex/250mg",
            license="CC0-1.0",
        ),
        Offer(
            inchikey=acid_nm.inchikey,
            vendor="ChemCo",
            catalog_id="VIAL-1G",
            pack_size=Quantity(value=1.0, unit="g", raw_unit_str="1 g"),
            price=45.0,
            currency="USD",
            price_basis="listed_catalog",
            source_id="local_stock",
            source_url="https://ex/1g",
            license="CC0-1.0",
        ),
        Offer(
            inchikey=acid_nm.inchikey,
            vendor="ChemCo",
            catalog_id="BOTTLE-5G",
            pack_size=Quantity(value=5.0, unit="g", raw_unit_str="5 g"),
            price=120.0,
            currency="USD",
            price_basis="listed_catalog",
            source_id="local_stock",
            source_url="https://ex/5g",
            license="CC0-1.0",
        ),
    ]

    # 1. Mixed-SKU optimum for 1.20 g: 1x 1g ($45) + 1x 250mg ($18) = $63.00 < single-SKU $90 (2x1g or 5x250mg)
    min_cost, plan = optimize_discrete_pack_order(offers, required_grams=1.20, currency="USD")
    assert min_cost == pytest.approx(63.0)
    assert sum(eff_g * cnt for _, eff_g, cnt in plan) >= 1.20 - 1e-9

    # 2. Shared-synthon pooling across two branches of a symmetrical bis-amide (N,N'-(1,4-phenylene)dibenzamide)
    target = "O=C(Nc1ccc(NC(=O)c2ccccc2)cc1)c1ccccc1"
    diamine = "Nc1ccc(N)cc1"
    inter = "Nc1ccc(NC(=O)c2ccccc2)cc1"
    diamine_nm = normalize_smiles(diamine)
    assert diamine_nm is not None

    # Target MW = 316.36, Benzoic acid MW = 122.12 (~0.386 g per branch -> 2 branches = 0.772 g total!)
    # Unpooled: each branch needs 0.386 g -> buys 2x 250mg ($36) per branch = $72 total (or $45 x 2 = $90 if only 1g vials)!
    # Pooled: 0.772 g total -> buys 1x 1g vial = $45 total -> saves $27!
    steps = [
        RouteStep(product=target, reactants=[acid, inter], template_or_reference="s1", source="t"),
        RouteStep(product=inter, reactants=[acid, diamine], template_or_reference="s2", source="t"),
    ]
    leaf_offers = {
        acid_nm.inchikey: offers,
        diamine_nm.inchikey: [
            Offer(
                inchikey=diamine_nm.inchikey,
                vendor="ChemCo",
                catalog_id="DIA-1G",
                pack_size=Quantity(value=1.0, unit="g", raw_unit_str="1 g"),
                price=10.0,
                currency="USD",
                price_basis="listed_catalog",
                source_id="local_stock",
                source_url="https://ex/dia",
                license="CC0-1.0",
            )
        ],
    }
    moq_res = optimize_route_moq_outlay(
        target,
        [acid, acid, diamine],
        leaf_offers,
        steps=steps,
        target_scale_g=1.0,
    )
    assert moq_res.complete is True
    assert moq_res.shared_leaf_count == 1
    assert moq_res.pooled_outlay is not None and moq_res.unpooled_outlay is not None
    assert moq_res.pooled_outlay.value == pytest.approx(55.0)     # $45 (1x 1g acid) + $10 (1x 1g diamine)
    assert moq_res.unpooled_outlay.value == pytest.approx(82.0)   # $36 + $36 (2x250mg per branch) + $10
    assert moq_res.pooling_savings_usd == pytest.approx(27.0)
