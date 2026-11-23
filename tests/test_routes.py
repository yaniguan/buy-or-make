"""Tests for retrosynthesis routes: acyclic `!@` SMARTS, atom conservation, spectators, 2-step lookup, and AiZynth."""

from __future__ import annotations

import json
import sqlite3
from collections import Counter
from pathlib import Path

import pytest
from rdkit import Chem

from buyormake import cost_aware_rank
from buyormake.canonical import normalize_smiles
from buyormake.routes import (
    RETRO_TEMPLATE_DELTAS,
    AiZynthFinderBackend,
    TemplateFallbackBackend,
    _walk_aizynth_tree,
    find_routes,
    parse_reaction_smiles_for_indexing,
)
from buyormake.schema import Offer, Route, RouteStep


def _elem_counts(smi: str) -> Counter[str]:
    m = Chem.MolFromSmiles(smi)
    assert m is not None
    return Counter(a.GetSymbol() for a in m.GetAtoms() if a.GetAtomicNum() > 1)


def test_acyclic_smarts_prevents_ring_cleavage_atom_duplication() -> None:
    tb = TemplateFallbackBackend()

    # 1. gamma-butyrolactone (C4H6O2): must NEVER split into bimolecular CCC(=O)O + CCO (5 carbons!)
    gbl_steps = tb.one_step("O=C1CCCO1")
    for st in gbl_steps:
        assert len(st.reactants) == 1
        assert sum(_elem_counts(r)["C"] for r in st.reactants) == 4
    assert any(st.reactants == ["O=C(O)CCCO"] for st in gbl_steps)

    # 2. Analog-6 (spiro-pyrimidine lactam, C12H16N4O): must NEVER split the 5-membered lactam into
    #    13-carbon CCC(=O)O + CC1(N)CCN(c2ncccn2)CC1; bimolecular disconnection must conserve C12!
    analog6 = "O=C1NC2(CC1)CCN(c1ncccn1)CC2"
    a6_steps = tb.one_step(analog6)
    assert len(a6_steps) >= 1
    for st in a6_steps:
        total_c = sum(_elem_counts(r)["C"] for r in st.reactants)
        total_n = sum(_elem_counts(r)["N"] for r in st.reactants)
        assert total_c == 12, f"Carbon duplication in {st.template_or_reference}: {st.reactants}"
        assert total_n == 4
    bimol_sets = [set(st.reactants) for st in a6_steps if len(st.reactants) == 2]
    assert {"Brc1ncccn1", "O=C1CCC2(CCNCC2)N1"} in bimol_sets

    # 3. Cyclic sulfonamide, cyclic urea, cyclic carbamate, THF, piperidine: zero atom-duplicating splits
    for cyc in ("O=S1(=O)CCCN1", "O=C1NCCN1", "O=C1OCCN1", "C1CCOC1", "C1CCNCC1"):
        p_counts = _elem_counts(cyc)
        for st in tb.one_step(cyc):
            tpl_name = st.template_or_reference.split(":", 1)[0]
            expected = p_counts + Counter(RETRO_TEMPLATE_DELTAS[tpl_name])
            actual = sum((_elem_counts(r) for r in st.reactants), Counter())
            assert actual == expected


def test_reaction_lookup_backend_and_multistep_expansion(tmp_path: Path) -> None:
    db_path = tmp_path / "reactions.sqlite"
    conn = sqlite3.connect(str(db_path))
    conn.execute(
        "CREATE TABLE reactions (id INTEGER PRIMARY KEY, product_inchikey TEXT, "
        "product_smiles TEXT, reactants_json TEXT, reaction_smiles TEXT, reference TEXT, source TEXT)"
    )
    # 2-step synthetic sequence using cubane/adamantyl-style bonds that TemplateFallback will not solve:
    # Step 1 (ORD): target -> intermediate + leaf_a
    # Step 2 (USPTO): intermediate -> leaf_b + leaf_c
    target = "CC12CC3CC(C1)(C3)C2C(=O)C1CC1"  # Analog-5 cage ketone
    intermediate = "CC12CC3CC(C1)(C3)C2C#N"
    leaf_a = "BrC1CC1"
    leaf_b = "CC12CC3CC(C1)(C3)C2Br"
    leaf_c = "[C-]#N.[Na+]"

    t_nm = normalize_smiles(target)
    i_nm = normalize_smiles(intermediate)
    assert t_nm and i_nm

    conn.execute(
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, reaction_smiles, reference, source) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        (
            t_nm.inchikey,
            t_nm.canonical_smiles,
            json.dumps([intermediate, leaf_a]),
            f"{intermediate}.{leaf_a}>>{target}",
            "ord-step1",
            "ord_lookup",
        ),
    )
    conn.execute(
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, reaction_smiles, reference, source) "
        "VALUES (?, ?, ?, ?, ?, ?)",
        (
            i_nm.inchikey,
            i_nm.canonical_smiles,
            json.dumps([leaf_b, leaf_c]),
            f"{leaf_b}.{leaf_c}>>{intermediate}",
            "US9999999",
            "uspto_lookup",
        ),
    )
    conn.commit()
    conn.close()

    stock = tmp_path / "stock.smi"
    stock.write_text(f"{leaf_a} cpr_br\n{leaf_b} cage_br\n{leaf_c} nacn\n", encoding="utf-8")

    # At max_steps=1, intermediate is not in stock -> unknown
    res_1 = cost_aware_rank(
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
    assert res_1[0].status == "unknown"

    # At max_steps=2, ReactionLookupBackend expands intermediate -> makeable in 2 steps!
    res_2 = cost_aware_rank(
        [target],
        max_steps=2,
        backends=("local_stock",),
        route_backend="lookup",
        offline=True,
        cache_path=tmp_path / "cache2.sqlite",
        stock_paths=[stock],
        reaction_db_path=db_path,
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    assert res_2[0].status == "makeable"
    assert res_2[0].n_steps == 2


def test_mapped_suzuki_spectator_filtering_in_indexing() -> None:
    # Mapped Suzuki coupling with unmapped THF (C1CCOC1), DMF (CN(C)C=O), Cs2CO3, and Pd(PPh3)4
    rxn_smi = (
        "[Br:99][c:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1."
        "[OH:97][B:98]([OH:96])[c:7]1[n:8][cH:9][cH:10][cH:11][cH:12]1."
        "C1CCOC1.CN(C)C=O.[Cs+].[Cs+].[O-]C(=O)[O-].[Pd](P(c1ccccc1)(c1ccccc1)c1ccccc1)"
        "(P(c1ccccc1)(c1ccccc1)c1ccccc1)(P(c1ccccc1)(c1ccccc1)c1ccccc1)P(c1ccccc1)(c1ccccc1)c1ccccc1"
        ">>[c:1]1([c:7]2[n:8][cH:9][cH:10][cH:11][cH:12]2)[cH:2][cH:3][cH:4][cH:5][cH:6]1"
    )
    extracted = parse_reaction_smiles_for_indexing(rxn_smi)
    assert len(extracted) == 1
    prod_nm, stoich_reactants, agents = extracted[0]
    assert prod_nm.canonical_smiles == "c1ccc(-c2ccccn2)cc1"
    # Only bromobenzene and pyridin-2-ylboronic acid are stoichiometric reactants!
    assert set(stoich_reactants) == {"Brc1ccccc1", "OB(O)c1ccccn1"}
    # THF, DMF, Cs2CO3, and Pd(PPh3)4 are classified as spectator agents!
    assert "C1CCOC1" in agents
    assert "CN(C)C=O" in agents


def test_aizynth_score_preservation_multi_or_node_and_auto_fallthrough() -> None:
    # 1. Multi-reaction OR-node tree: branch 0 has out-of-stock leaf, branch 1 has in-stock leaves
    tree = {
        "type": "mol",
        "smiles": "O=C(Nc1ccccc1)c1ccccc1",
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "metadata": {"template_code": "bad_branch"},
                "children": [
                    {"type": "mol", "smiles": "O=C(Cl)c1ccccc1", "in_stock": False, "children": []},
                    {"type": "mol", "smiles": "Nc1ccccc1", "in_stock": True, "children": []},
                ],
            },
            {
                "type": "reaction",
                "metadata": {"template_code": "good_branch"},
                "children": [
                    {"type": "mol", "smiles": "O=C(O)c1ccccc1", "in_stock": True, "children": []},
                    {"type": "mol", "smiles": "Nc1ccccc1", "in_stock": True, "children": []},
                ],
            },
        ],
    }
    steps: list[RouteStep] = []
    leaves: list[str] = []
    leaves_ok: list[bool] = []
    _walk_aizynth_tree(tree, steps, leaves, leaves_ok)
    assert len(steps) == 1
    assert steps[0].template_or_reference == "good_branch"
    assert leaves_ok == [True, True]

    # 2. Score preservation in find_routes AND fallthrough in `auto` mode when AiZynth is unsolved
    target = "O=C(Nc1ccccc1)c1ccccc1"
    t_nm = normalize_smiles(target)
    acid_nm = normalize_smiles("O=C(O)c1ccccc1")
    amine_nm = normalize_smiles("Nc1ccccc1")
    assert t_nm and acid_nm and amine_nm

    stock_by_ik = {
        acid_nm.inchikey: [
            Offer(
                inchikey=acid_nm.inchikey,
                vendor="V",
                source_id="local_stock",
                source_url="https://example.org/a",
                license="CC0-1.0",
            )
        ],
        amine_nm.inchikey: [
            Offer(
                inchikey=amine_nm.inchikey,
                vendor="V",
                source_id="local_stock",
                source_url="https://example.org/b",
                license="CC0-1.0",
            )
        ],
    }

    class DummySolvedAZ(AiZynthFinderBackend):
        def is_ready(self) -> bool:
            return True

        def search(self, target_smiles: str, target_inchikey: str, max_steps: int) -> list[Route]:
            return [
                Route(
                    target_inchikey=target_inchikey,
                    n_steps=1,
                    steps=steps,
                    leaves=leaves,
                    leaves_purchasable=[True, True],
                    method="aizynthfinder (MCTS retrosynthesis)",
                    score=0.9375,
                    est_cost=None,
                    est_cost_method="leaf sum",
                    est_cost_complete=False,
                )
            ]

    solved_routes = find_routes(
        t_nm.canonical_smiles,
        t_nm.inchikey,
        max_steps=2,
        route_backend="aizynthfinder",
        offer_lookup=lambda ik, smi: stock_by_ik.get(ik, []),
        aizynth_backend=DummySolvedAZ(),
    )
    assert len(solved_routes) == 1
    assert solved_routes[0].score == pytest.approx(0.9375)
    assert all(solved_routes[0].leaves_purchasable)

    class DummyUnsolvedAZ(AiZynthFinderBackend):
        def is_ready(self) -> bool:
            return True

        def search(self, target_smiles: str, target_inchikey: str, max_steps: int) -> list[Route]:
            # Returns only an unsolved route (acid chloride not in stock_by_ik)
            return [
                Route(
                    target_inchikey=target_inchikey,
                    n_steps=1,
                    steps=[
                        RouteStep(
                            product=target_smiles,
                            reactants=["O=C(Cl)c1ccccc1", "Nc1ccccc1"],
                            template_or_reference="az_unsolved",
                            source="aizynthfinder",
                        )
                    ],
                    leaves=["O=C(Cl)c1ccccc1", "Nc1ccccc1"],
                    leaves_purchasable=[False, True],
                    method="aizynthfinder (MCTS retrosynthesis)",
                    score=0.42,
                    est_cost=None,
                    est_cost_method="leaf sum",
                    est_cost_complete=False,
                )
            ]

    # In `auto` mode, unsolved AiZynth route MUST NOT block fallback to TemplateFallbackBackend!
    auto_routes = find_routes(
        t_nm.canonical_smiles,
        t_nm.inchikey,
        max_steps=2,
        route_backend="auto",
        offer_lookup=lambda ik, smi: stock_by_ik.get(ik, []),
        aizynth_backend=DummyUnsolvedAZ(),
    )
    assert any(r.leaves and all(r.leaves_purchasable) for r in auto_routes)
    assert auto_routes[0].steps[0].source == "template_fallback"
