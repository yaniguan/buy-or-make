#!/usr/bin/env python3
"""Reproducible leakage-free CPU multi-step & single-step retrosynthesis benchmark on public corpora.

Evaluates `buy-or-make` across:
1. `PaRoutes 2.0` (`n1` and `n5` stopping stocks + multi-step patent route targets, Zenodo `7341155`, `CC-BY-4.0`)
2. `SynPlanner` SAscore-stratified retrosynthesis benchmark targets + building blocks (`MIT`)
3. `Syntheseus USPTO-50K` held-out `raw_test.csv` single-step top-k accuracy (`MIT`, train+val only indexed)

Strict Split & Target-Exclusion Hygiene (Issue #35):
- Every reaction extracted from `PaRoutes 2.0` (`ref_routes_n1.json`, `ref_routes_n5.json`) is tagged with
  its `tree_index`.
- For all leakage-free multi-step methods (`lookup_1step_smarts_2step`, `retro_star_max4_leakage_free`,
  headline `retro_star_max5_leakage_free`, and `retro_star_max5_uspto50k_only`), the reaction lookup SQLite
  index strictly excludes:
  (a) ALL reactions from any of the `160` sampled `PaRoutes` evaluation trees (`80` `n1` + `80` `n5`), AND
  (b) ANY reaction across the entire index (held-out non-target `PaRoutes` trees AND `USPTO-50K` train+val)
      whose product `InChIKey` (27-char) or 14-char skeletal connectivity block matches any evaluated
      `PaRoutes` or `SynPlanner` target!
- A separate `retro_star_max5_retrieval_upper_bound` row retains all `PaRoutes` trees (including each
  target's own reference route) solely as an explicitly labeled retrieval upper bound, NEVER as the headline.

Writes `results/retrosynthesis_benchmark.json` and optionally syncs the benchmark section in `README.md`.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import platform
import random
import sqlite3
import statistics
import sys
import time
import zipfile
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

import rdkit

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from buyormake.canonical import file_sha256, normalize_smiles  # noqa: E402
from buyormake.open_corpora import (  # noqa: E402
    CORPUS_ARTIFACTS,
    DEFAULT_OPEN_CORPORA_DIR,
    _index_reaction_pairs,
    _init_reaction_db,
    download_corpus_artifact,
    import_paroutes_stock,
    import_synplanner_building_blocks,
    load_synplanner_targets,
    parse_paroutes_route_tree,
)
from buyormake.retro_star import EXTENDED_RETRO_STAR_TEMPLATES, RetroStarPlanner  # noqa: E402
from buyormake.routes import (  # noqa: E402
    RETRO_TEMPLATES,
    ReactionLookupBackend,
    TemplateFallbackBackend,
    _build_route,
    _is_purchasable,
    _score_route_for_selection,
    find_routes,
    parse_reaction_smiles_for_indexing,
)
from buyormake.schema import Offer, Route  # noqa: E402
from buyormake.scores import DEFAULT_SCSCORE_WEIGHTS, compute_sascore  # noqa: E402
from buyormake.sources import SOURCE_BY_ID  # noqa: E402
from buyormake.stock import LocalStockBackend  # noqa: E402

DEFAULT_RESULTS_JSON = REPO_ROOT / "results" / "retrosynthesis_benchmark.json"
README_PATH = REPO_ROOT / "README.md"
BENCHMARK_START_MARKER = "<!-- BENCHMARK_TABLE_START -->"
BENCHMARK_END_MARKER = "<!-- BENCHMARK_TABLE_END -->"

HEADLINE_METHOD = "retro_star_max5_leakage_free"
N_CORE_TEMPLATES = len(RETRO_TEMPLATES)
N_EXT_TEMPLATES = len(EXTENDED_RETRO_STAR_TEMPLATES)

REQUIRED_ARTIFACTS = (
    "paroutes_stock_n1",
    "paroutes_stock_n5",
    "paroutes_routes_n1",
    "paroutes_routes_n5",
    "syntheseus_uspto_50k",
    "synplanner_building_blocks",
    "synplanner_targets_sa_1_5_2_5",
    "synplanner_targets_sa_2_5_3_5",
    "synplanner_targets_sa_3_5_4_5",
    "synplanner_targets_sa_4_5_5_5",
)


def _ensure_artifacts(corpora_dir: Path) -> dict[str, dict[str, Any]]:
    """Download/verify all required public benchmark files and compute exact SHA-256 digests."""
    manifest: dict[str, dict[str, Any]] = {}
    for key in REQUIRED_ARTIFACTS:
        spec = CORPUS_ARTIFACTS[key]
        src = SOURCE_BY_ID[spec.source_id]
        path = download_corpus_artifact(key, out_root=corpora_dir)
        digest = file_sha256(path)
        manifest[key] = {
            "artifact_key": key,
            "source_id": spec.source_id,
            "filename": path.name,
            "spdx_license": src.spdx_license or src.license,
            "url": spec.url,
            "size_bytes": path.stat().st_size,
            "sha256": digest,
        }

    sc_path = REPO_ROOT / DEFAULT_SCSCORE_WEIGHTS
    if not sc_path.exists():
        import scripts.fetch_scscore as fetch_scscore

        fetch_scscore.main(["--out", str(sc_path)])
    manifest["scscore_weights"] = {
        "artifact_key": "scscore_weights",
        "source_id": "scscore",
        "filename": sc_path.name,
        "spdx_license": "MIT",
        "url": "https://raw.githubusercontent.com/connorcoley/scscore/master/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz",
        "size_bytes": sc_path.stat().st_size,
        "sha256": file_sha256(sc_path),
    }
    return manifest


def _is_convergent_route(route: Route) -> bool:
    prod_set = {s.product for s in route.steps}
    return any(sum(1 for r in s.reactants if r in prod_set) >= 2 for s in route.steps)


def _sa_band_label(sa: float) -> str:
    if sa < 2.5:
        return "[1.5, 2.5)"
    if sa < 3.5:
        return "[2.5, 3.5)"
    if sa < 4.5:
        return "[3.5, 4.5)"
    return "[4.5, 10.0]"


def _sample_paroutes_targets(
    routes_json_path: Path,
    *,
    dataset_label: str,
    per_depth: int,
    depths: tuple[int, ...] = (2, 3, 4, 5),
    seed: int = 20261127,
) -> tuple[list[dict[str, Any]], list[tuple[int, str, str]], set[int], int]:
    """Deterministically sample stratified patent targets by route depth `{2,3,4,5}` and return
    tree-index-tagged `(tree_index, rxn_smiles, ref)` tuples plus `sampled_tree_indices`.
    """
    raw_trees = json.loads(routes_json_path.read_text(encoding="utf-8"))
    total_trees = len(raw_trees)
    by_depth: dict[int, list[dict[str, Any]]] = {d: [] for d in depths}
    seen_target_iks: set[str] = set()
    sampled_tree_indices: set[int] = set()

    # Deterministically shuffle indices to avoid file-ordering bias before parsing stratified targets
    rng = random.Random(seed + (101 if "n1" in dataset_label else 202))
    indices = list(range(total_trees))
    rng.shuffle(indices)

    for idx in indices:
        if all(len(by_depth[d]) >= per_depth for d in depths):
            break
        route = parse_paroutes_route_tree(raw_trees[idx])
        if route is None or route.n_steps not in by_depth:
            continue
        if len(by_depth[route.n_steps]) >= per_depth:
            continue
        if route.target_inchikey in seen_target_iks:
            continue
        if not (route.leaves and all(route.leaves_purchasable)):
            continue
        root_smi = route.steps[0].product
        nm = normalize_smiles(root_smi)
        if nm is None:
            continue
        sa_val = compute_sascore(nm.canonical_smiles, mol=nm.mol) or 2.5
        seen_target_iks.add(nm.inchikey)
        sampled_tree_indices.add(idx)
        by_depth[route.n_steps].append(
            {
                "smiles": nm.canonical_smiles,
                "inchikey": nm.inchikey,
                "dataset": dataset_label,
                "tree_index": idx,
                "ref_depth": route.n_steps,
                "sascore": round(float(sa_val), 4),
                "sa_band": _sa_band_label(sa_val),
                "is_convergent_ref": _is_convergent_route(route),
            }
        )

    tagged_rxns: list[tuple[int, str, str]] = []
    seen_tagged: set[tuple[int, str, str]] = set()

    def _collect_tree_reactions(node: dict[str, Any], tree_idx: int, prefix: str) -> None:
        ntype = node.get("type", "mol")
        if ntype == "mol":
            p_raw = str(node.get("smiles") or "").strip()
            for ch in node.get("children") or []:
                if isinstance(ch, dict) and ch.get("type") == "reaction":
                    mol_ch = [
                        mc
                        for mc in (ch.get("children") or [])
                        if isinstance(mc, dict) and mc.get("type", "mol") == "mol"
                    ]
                    r_smis = [str(mc.get("smiles") or "").strip() for mc in mol_ch if mc.get("smiles")]
                    if p_raw and r_smis:
                        meta = ch.get("metadata") if isinstance(ch.get("metadata"), dict) else {}
                        ref = str(meta.get("ID") or meta.get("id") or prefix)
                        rsmi = f"{'.'.join(r_smis)}>>{p_raw}"
                        sig = (tree_idx, rsmi, ref)
                        if sig not in seen_tagged:
                            seen_tagged.add(sig)
                            tagged_rxns.append(sig)
                    for mc in mol_ch:
                        _collect_tree_reactions(mc, tree_idx, prefix)

    for idx, item in enumerate(raw_trees):
        root = item.get("dict") if isinstance(item.get("dict"), dict) else item
        if isinstance(root, dict):
            _collect_tree_reactions(root, idx, f"{dataset_label}:r{idx}")

    targets: list[dict[str, Any]] = []
    for d in depths:
        targets.extend(by_depth[d])
    return targets, tagged_rxns, sampled_tree_indices, total_trees


def _sample_synplanner_targets(
    corpora_dir: Path,
    *,
    per_band: int,
    seed: int = 20261127,
) -> list[dict[str, Any]]:
    band_specs = (
        ("[1.5, 2.5)", "targets_with_sascore_1.5_2.5.smi"),
        ("[2.5, 3.5)", "targets_with_sascore_2.5_3.5.smi"),
        ("[3.5, 4.5)", "targets_with_sascore_3.5_4.5.smi"),
        ("[4.5, 10.0]", "targets_with_sascore_4.5_5.5.smi"),
    )
    rng = random.Random(seed + 303)
    out: list[dict[str, Any]] = []
    seen_iks: set[str] = set()
    for band_label, fname in band_specs:
        p = corpora_dir / "synplanner_data" / fname
        all_smis = load_synplanner_targets(p)
        shuffled = list(all_smis)
        rng.shuffle(shuffled)
        picked = 0
        for smi in shuffled:
            nm = normalize_smiles(smi)
            if nm is None or nm.inchikey in seen_iks:
                continue
            seen_iks.add(nm.inchikey)
            sa_val = compute_sascore(nm.canonical_smiles, mol=nm.mol) or 3.0
            out.append(
                {
                    "smiles": nm.canonical_smiles,
                    "inchikey": nm.inchikey,
                    "dataset": "synplanner_sascore",
                    "tree_index": None,
                    "ref_depth": None,
                    "sascore": round(float(sa_val), 4),
                    "sa_band": band_label,
                    "is_convergent_ref": None,
                }
            )
            picked += 1
            if picked >= per_band:
                break
    return out


def _index_filtered_reaction_pairs(
    pairs: Iterable[tuple[Optional[str], Optional[str]]],
    *,
    db_path: Path,
    step_source: str,
    default_prefix: str,
    excluded_product_inchikeys: frozenset[str] = frozenset(),
    excluded_connectivity_blocks: frozenset[str] = frozenset(),
) -> tuple[int, int]:
    """Index reactions into `db_path`, skipping any reaction whose product matches an excluded target.

    Returns `(inserted_rows, skipped_target_product_rows)`.
    """
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path), timeout=30.0)
    _init_reaction_db(conn)
    sql = (
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
        "agents_json, reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?, ?)"
    )
    batch: list[tuple[str, str, str, str, str, str, str]] = []
    inserted = 0
    skipped_target_prod = 0
    try:
        for idx, (raw_rxn, ref) in enumerate(pairs):
            if not raw_rxn:
                continue
            rxn_smi = str(raw_rxn).strip()
            if rxn_smi.count(">") == 2 and ">>" not in rxn_smi:
                lhs, mid, rhs = rxn_smi.split(">")
                if mid.strip():
                    rxn_for_parse = f"{lhs.strip()}>{mid.strip()}>>{rhs.strip()}"
                else:
                    rxn_for_parse = f"{lhs.strip()}>>{rhs.strip()}"
            else:
                rxn_for_parse = rxn_smi
            if ">>" not in rxn_for_parse:
                continue
            ref_str = str(ref).strip() if ref and str(ref).strip() else f"{default_prefix}:{idx}"
            for p_nm, stoich_reactants, agents in parse_reaction_smiles_for_indexing(rxn_for_parse):
                if (
                    p_nm.inchikey in excluded_product_inchikeys
                    or p_nm.connectivity_block in excluded_connectivity_blocks
                ):
                    skipped_target_prod += 1
                    continue
                batch.append(
                    (
                        p_nm.inchikey,
                        p_nm.canonical_smiles,
                        json.dumps(stoich_reactants),
                        json.dumps(agents),
                        rxn_smi,
                        ref_str,
                        step_source,
                    )
                )
                inserted += 1
                if len(batch) >= 2000:
                    conn.executemany(sql, batch)
                    conn.commit()
                    batch.clear()
        if batch:
            conn.executemany(sql, batch)
            conn.commit()
    finally:
        conn.close()
    return inserted, skipped_target_prod


def verify_reaction_index_hygiene(
    db_path: Path,
    *,
    paroutes_target_inchikeys: set[str],
    all_target_inchikeys: set[str],
) -> dict[str, int]:
    """Query `db_path` directly to measure exact InChIKey & 14-char connectivity overlap with targets."""
    conn = sqlite3.connect(str(db_path))
    try:
        total_rows = int(conn.execute("SELECT COUNT(*) FROM reactions").fetchone()[0])
        distinct_iks = {
            str(r[0]) for r in conn.execute("SELECT DISTINCT product_inchikey FROM reactions").fetchall()
        }
    finally:
        conn.close()

    paroutes_overlap = len(distinct_iks & paroutes_target_inchikeys)
    all_target_overlap = len(distinct_iks & all_target_inchikeys)
    conn_blocks_db = {ik[:14] for ik in distinct_iks if len(ik) >= 14}
    conn_blocks_targets = {ik[:14] for ik in all_target_inchikeys if len(ik) >= 14}
    connectivity_overlap = len(conn_blocks_db & conn_blocks_targets)
    return {
        "indexed_reaction_rows": total_rows,
        "distinct_product_inchikeys": len(distinct_iks),
        "paroutes_target_product_inchikey_overlap": paroutes_overlap,
        "all_target_product_inchikey_overlap": all_target_overlap,
        "all_target_connectivity_block_overlap": connectivity_overlap,
    }


def _build_benchmark_databases(
    corpora_dir: Path,
    cache_dir: Path,
    *,
    tagged_rxn_n1: list[tuple[int, str, str]],
    sampled_trees_n1: set[int],
    tagged_rxn_n5: list[tuple[int, str, str]],
    sampled_trees_n5: set[int],
    evaluated_target_inchikeys: set[str],
) -> tuple[dict[str, Path], dict[str, int]]:
    """Build stopping-stock DBs and strictly leakage-free vs retrieval-upper-bound reaction indexes."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    stock_n1_db = cache_dir / "bench_stock_n1.sqlite"
    stock_n5_db = cache_dir / "bench_stock_n5.sqlite"
    stock_sp_db = cache_dir / "bench_stock_synplanner.sqlite"

    # 1. Single-step USPTO-50K train+val DB (raw_test.csv strictly excluded)
    rxn_uspto_single_step_db = cache_dir / "bench_rxn_uspto50k_trainval_singlestep.sqlite"
    # 2. Multi-step strict zero-PaRoutes DB (USPTO-50K train+val minus any evaluated target products)
    rxn_uspto_multistep_db = cache_dir / "bench_rxn_uspto50k_strict_leakfree.sqlite"
    # 3. Headline leakage-free multi-step DB (~19,840 non-target PaRoutes trees + USPTO-50K train+val, 0 target products)
    rxn_heldout_leakfree_db = cache_dir / "bench_rxn_heldout_paroutes_uspto50k_leakfree.sqlite"
    # 4. Reference-route retrieval upper-bound DB (includes sampled targets' own PaRoutes trees)
    rxn_retrieval_upper_db = cache_dir / "bench_rxn_retrieval_upper_bound.sqlite"

    if not stock_n1_db.exists():
        import_paroutes_stock(
            corpora_dir / "paroutes" / "stock_n1.txt", db_path=stock_n1_db, subset="n1"
        )
    if not stock_n5_db.exists():
        import_paroutes_stock(
            corpora_dir / "paroutes" / "stock_n5.txt", db_path=stock_n5_db, subset="n5"
        )
    if not stock_sp_db.exists():
        import_synplanner_building_blocks(
            corpora_dir / "synplanner_data" / "building_blocks_em_sa_ln.smi",
            db_path=stock_sp_db,
        )
        import_paroutes_stock(
            corpora_dir / "paroutes" / "stock_n1.txt", db_path=stock_sp_db, subset="n1"
        )

    excluded_iks = frozenset(evaluated_target_inchikeys)
    excluded_blocks = frozenset(ik[:14] for ik in evaluated_target_inchikeys if len(ik) >= 14)

    uspto_zip = corpora_dir / "syntheseus_uspto_50k" / "uspto50k.zip"

    def _iter_uspto_trainval() -> Iterable[tuple[Optional[str], Optional[str]]]:
        with zipfile.ZipFile(uspto_zip, "r") as zf:
            for split_csv in ("raw_train.csv", "raw_val.csv"):
                with zf.open(split_csv, "r") as raw:
                    rdr = csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"))
                    for i, row in enumerate(rdr):
                        rxn = row.get("reactants>reagents>production") or row.get("reactions")
                        ref = row.get("id") or f"{split_csv}:{i}"
                        yield rxn, ref

    if not rxn_uspto_single_step_db.exists():
        _index_reaction_pairs(
            _iter_uspto_trainval(),
            db_path=rxn_uspto_single_step_db,
            step_source="syntheseus_uspto_50k_lookup",
            default_prefix="uspto50k",
        )

    # Always rebuild target-filtered multi-step reaction indexes if missing so target exclusions are exact
    filter_counts = {
        "paroutes_target_tree_reactions_excluded": 0,
        "heldout_paroutes_target_product_reactions_excluded": 0,
        "uspto50k_target_product_reactions_excluded": 0,
    }

    heldout_paroutes_pairs: list[tuple[str, str]] = []
    seen_heldout_pairs: set[tuple[str, str]] = set()
    all_paroutes_pairs: list[tuple[str, str]] = []
    seen_all_pairs: set[tuple[str, str]] = set()

    for tagged_list, sampled_set in (
        (tagged_rxn_n1, sampled_trees_n1),
        (tagged_rxn_n5, sampled_trees_n5),
    ):
        for tree_idx, rsmi, ref in tagged_list:
            pair = (rsmi, ref)
            if pair not in seen_all_pairs:
                seen_all_pairs.add(pair)
                all_paroutes_pairs.append(pair)
            if tree_idx in sampled_set:
                filter_counts["paroutes_target_tree_reactions_excluded"] += 1
            else:
                if pair not in seen_heldout_pairs:
                    seen_heldout_pairs.add(pair)
                    heldout_paroutes_pairs.append(pair)

    if rxn_uspto_multistep_db.exists():
        rxn_uspto_multistep_db.unlink()
    _, skipped_u = _index_filtered_reaction_pairs(
        _iter_uspto_trainval(),
        db_path=rxn_uspto_multistep_db,
        step_source="syntheseus_uspto_50k_lookup",
        default_prefix="uspto50k",
        excluded_product_inchikeys=excluded_iks,
        excluded_connectivity_blocks=excluded_blocks,
    )
    filter_counts["uspto50k_target_product_reactions_excluded"] = skipped_u

    if rxn_heldout_leakfree_db.exists():
        rxn_heldout_leakfree_db.unlink()
    conn_src = sqlite3.connect(str(rxn_uspto_multistep_db))
    conn_dst = sqlite3.connect(str(rxn_heldout_leakfree_db))
    conn_src.backup(conn_dst)
    conn_dst.close()
    conn_src.close()
    _, skipped_p = _index_filtered_reaction_pairs(
        heldout_paroutes_pairs,
        db_path=rxn_heldout_leakfree_db,
        step_source="paroutes_lookup",
        default_prefix="paroutes_heldout",
        excluded_product_inchikeys=excluded_iks,
        excluded_connectivity_blocks=excluded_blocks,
    )
    filter_counts["heldout_paroutes_target_product_reactions_excluded"] = skipped_p

    if not rxn_retrieval_upper_db.exists():
        conn_src_u = sqlite3.connect(str(rxn_uspto_single_step_db))
        conn_dst_u = sqlite3.connect(str(rxn_retrieval_upper_db))
        conn_src_u.backup(conn_dst_u)
        conn_dst_u.close()
        conn_src_u.close()
        _index_reaction_pairs(
            all_paroutes_pairs,
            db_path=rxn_retrieval_upper_db,
            step_source="paroutes_lookup",
            default_prefix="paroutes",
        )

    return (
        {
            "stock_n1": stock_n1_db,
            "stock_n5": stock_n5_db,
            "stock_synplanner": stock_sp_db,
            "rxn_uspto_single_step": rxn_uspto_single_step_db,
            "rxn_uspto_strict_leakfree": rxn_uspto_multistep_db,
            "rxn_heldout_leakfree": rxn_heldout_leakfree_db,
            "rxn_retrieval_upper_bound": rxn_retrieval_upper_db,
        },
        filter_counts,
    )


def _make_offer_lookup(stock_backend: LocalStockBackend) -> Callable[[str, str], list[Offer]]:
    memo: dict[str, list[Offer]] = {}

    def _lookup(ik: str, canonical_smiles: str) -> list[Offer]:
        if ik in memo:
            return memo[ik]
        offs = list(stock_backend.lookup(ik))
        if not offs:
            offs = stock_backend.lookup_parent_equivalents(ik, exclude_inchikey=ik)
        memo[ik] = offs
        return offs

    return _lookup


def _eval_baseline_1step_lookup_2step_smarts(
    smiles: str,
    inchikey: str,
    offer_lookup: Callable[[str, str], list[Offer]],
    lookup_backend: ReactionLookupBackend,
    template_backend: TemplateFallbackBackend,
) -> list[Route]:
    """Method 2: 1-step held-out reaction lookup + 2-step greedy SMARTS (`17` core templates)."""
    offer_cache: dict[str, list[Offer]] = {}
    routes: list[Route] = []

    for st1 in lookup_backend.one_step(inchikey, smiles, limit=10):
        r1 = _build_route(
            smiles,
            inchikey,
            [st1],
            list(st1.reactants),
            f"reaction_lookup_1step ({st1.source})",
            offer_cache,
            offer_lookup,
        )
        routes.append(r1)

    if any(r.leaves and all(r.leaves_purchasable) for r in routes):
        routes.sort(key=_score_route_for_selection)
        return routes

    t_routes = find_routes(
        smiles,
        inchikey,
        max_steps=2,
        route_backend="template",
        planner="lookup_smarts",
        offer_lookup=offer_lookup,
        template_backend=template_backend,
        offer_cache=offer_cache,
    )
    routes.extend(t_routes)
    routes.sort(key=_score_route_for_selection)
    return routes


def _summarize_method_runs(records: list[dict[str, Any]]) -> dict[str, Any]:
    n = len(records)
    if n == 0:
        return {
            "n_targets": 0,
            "solved_count": 0,
            "solve_rate": 0.0,
            "median_time_ms": 0.0,
            "mean_time_ms": 0.0,
            "mean_leaf_stock_completeness": 0.0,
            "solved_leaf_stock_completeness": 1.0,
            "mean_steps_when_solved": None,
        }
    solved = [r for r in records if r["solved"]]
    times = [float(r["time_ms"]) for r in records]
    leaf_comp = [float(r["leaf_completeness"]) for r in records]
    solved_steps = [int(r["n_steps"]) for r in solved if r["n_steps"] is not None]
    return {
        "n_targets": n,
        "solved_count": len(solved),
        "solve_rate": round(len(solved) / n, 4),
        "median_time_ms": round(float(statistics.median(times)), 3),
        "mean_time_ms": round(float(statistics.mean(times)), 3),
        "mean_leaf_stock_completeness": round(float(statistics.mean(leaf_comp)), 4),
        "solved_leaf_stock_completeness": (
            round(float(statistics.mean(r["leaf_completeness"] for r in solved)), 4)
            if solved
            else 1.0
        ),
        "mean_steps_when_solved": (
            round(float(statistics.mean(solved_steps)), 3) if solved_steps else None
        ),
    }


def _evaluate_uspto50k_single_step_topk(
    uspto_zip: Path,
    lookup_trainval: ReactionLookupBackend,
    *,
    n_test_samples: int = 500,
    seed: int = 20261127,
) -> dict[str, Any]:
    """Evaluate single-step top-k reactant-set recovery on held-out `USPTO-50K raw_test.csv`."""
    raw_rows: list[tuple[str, str]] = []
    with zipfile.ZipFile(uspto_zip, "r") as zf:
        with zf.open("raw_test.csv", "r") as raw:
            rdr = csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"))
            for row in rdr:
                rxn = row.get("reactants>reagents>production") or row.get("reactions")
                if rxn:
                    raw_rows.append((row.get("id") or "", rxn))

    rng = random.Random(seed + 505)
    indices = list(range(len(raw_rows)))
    rng.shuffle(indices)

    test_items: list[tuple[str, str, Any, frozenset[str]]] = []
    for idx in indices:
        if len(test_items) >= n_test_samples:
            break
        _, rxn_smi = raw_rows[idx]
        parsed = parse_reaction_smiles_for_indexing(rxn_smi)
        if not parsed:
            continue
        p_nm, gt_reactants, _ = parsed[0]
        if not gt_reactants or p_nm.mol is None:
            continue
        test_items.append((p_nm.canonical_smiles, p_nm.inchikey, p_nm.mol, frozenset(gt_reactants)))

    tb_17 = TemplateFallbackBackend(max_branch_per_node=30)
    rs_ext = RetroStarPlanner(max_branch_per_node=40)

    ks = (1, 3, 5, 10)
    base_hits = {k: 0 for k in ks}
    ext_hits = {k: 0 for k in ks}
    base_times: list[float] = []
    ext_times: list[float] = []

    for smi, ik, mol, gt_set in test_items:
        t0 = time.perf_counter()
        b_preds: list[frozenset[str]] = []
        seen_b: set[frozenset[str]] = set()
        for st in lookup_trainval.one_step(ik, smi, limit=10):
            fs = frozenset(st.reactants)
            if fs not in seen_b:
                seen_b.add(fs)
                b_preds.append(fs)
        for st in tb_17.one_step(smi, mol=mol):
            fs = frozenset(st.reactants)
            if fs not in seen_b:
                seen_b.add(fs)
                b_preds.append(fs)
        base_times.append((time.perf_counter() - t0) * 1000.0)
        for k in ks:
            if gt_set in b_preds[:k]:
                base_hits[k] += 1

        t1 = time.perf_counter()
        e_preds: list[frozenset[str]] = []
        seen_e: set[frozenset[str]] = set()
        for st in lookup_trainval.one_step(ik, smi, limit=10):
            fs = frozenset(st.reactants)
            if fs not in seen_e:
                seen_e.add(fs)
                e_preds.append(fs)
        for _, st in rs_ext.expand_extended_smarts(smi, mol=mol, allow_unstocked_fgi=True):
            fs = frozenset(st.reactants)
            if fs not in seen_e:
                seen_e.add(fs)
                e_preds.append(fs)
        ext_times.append((time.perf_counter() - t1) * 1000.0)
        for k in ks:
            if gt_set in e_preds[:k]:
                ext_hits[k] += 1

    total = len(test_items)
    ext_key = f"retro_star_lookup_trainval_plus_{N_EXT_TEMPLATES}_smarts"
    return {
        "dataset": "syntheseus_uspto_50k_test",
        "split_hygiene": "raw_train.csv + raw_val.csv indexed; raw_test.csv strictly held out (0 test leakage)",
        "n_test_targets": total,
        "n_extended_smarts_templates": N_EXT_TEMPLATES,
        "methods": {
            "lookup_trainval_plus_17_smarts": {
                "top_1": round(base_hits[1] / total, 4),
                "top_3": round(base_hits[3] / total, 4),
                "top_5": round(base_hits[5] / total, 4),
                "top_10": round(base_hits[10] / total, 4),
                "median_time_ms": round(float(statistics.median(base_times)), 3),
            },
            ext_key: {
                "top_1": round(ext_hits[1] / total, 4),
                "top_3": round(ext_hits[3] / total, 4),
                "top_5": round(ext_hits[5] / total, 4),
                "top_10": round(ext_hits[10] / total, 4),
                "median_time_ms": round(float(statistics.median(ext_times)), 3),
            },
        },
    }


def run_full_benchmark(
    *,
    corpora_dir: Path = DEFAULT_OPEN_CORPORA_DIR,
    cache_dir: Path = Path("data/benchmark_cache"),
    paroutes_per_depth: int = 20,
    synplanner_per_band: int = 20,
    uspto50k_test_n: int = 500,
    retro_star_expansions: int = 120,
    seed: int = 20261127,
) -> dict[str, Any]:
    manifest = _ensure_artifacts(corpora_dir)

    targets_n1, tagged_n1, sampled_trees_n1, total_trees_n1 = _sample_paroutes_targets(
        corpora_dir / "paroutes" / "ref_routes_n1.json",
        dataset_label="paroutes_n1",
        per_depth=paroutes_per_depth,
        seed=seed,
    )
    targets_n5, tagged_n5, sampled_trees_n5, total_trees_n5 = _sample_paroutes_targets(
        corpora_dir / "paroutes" / "ref_routes_n5.json",
        dataset_label="paroutes_n5",
        per_depth=paroutes_per_depth,
        seed=seed,
    )
    targets_sp = _sample_synplanner_targets(
        corpora_dir, per_band=synplanner_per_band, seed=seed
    )
    all_targets = [*targets_n1, *targets_n5, *targets_sp]
    assert len(all_targets) >= 200, f"Expected >= 200 multi-step targets, got {len(all_targets)}"

    paroutes_target_iks = {t["inchikey"] for t in [*targets_n1, *targets_n5]}
    all_target_iks = {t["inchikey"] for t in all_targets}

    db_paths, filter_counts = _build_benchmark_databases(
        corpora_dir,
        cache_dir,
        tagged_rxn_n1=tagged_n1,
        sampled_trees_n1=sampled_trees_n1,
        tagged_rxn_n5=tagged_n5,
        sampled_trees_n5=sampled_trees_n5,
        evaluated_target_inchikeys=all_target_iks,
    )

    # Hard runtime verification of zero target-route / target-product leakage on leakage-free DBs
    hygiene_heldout = verify_reaction_index_hygiene(
        db_paths["rxn_heldout_leakfree"],
        paroutes_target_inchikeys=paroutes_target_iks,
        all_target_inchikeys=all_target_iks,
    )
    hygiene_uspto = verify_reaction_index_hygiene(
        db_paths["rxn_uspto_strict_leakfree"],
        paroutes_target_inchikeys=paroutes_target_iks,
        all_target_inchikeys=all_target_iks,
    )
    hygiene_upper = verify_reaction_index_hygiene(
        db_paths["rxn_retrieval_upper_bound"],
        paroutes_target_inchikeys=paroutes_target_iks,
        all_target_inchikeys=all_target_iks,
    )
    assert hygiene_heldout["paroutes_target_product_inchikey_overlap"] == 0
    assert hygiene_heldout["all_target_product_inchikey_overlap"] == 0
    assert hygiene_heldout["all_target_connectivity_block_overlap"] == 0
    assert hygiene_uspto["all_target_product_inchikey_overlap"] == 0

    heldout_trees_indexed = (total_trees_n1 - len(sampled_trees_n1)) + (
        total_trees_n5 - len(sampled_trees_n5)
    )
    excluded_trees_total = len(sampled_trees_n1) + len(sampled_trees_n5)

    stock_n1 = LocalStockBackend([db_paths["stock_n1"]])
    stock_n5 = LocalStockBackend([db_paths["stock_n5"]])
    stock_sp = LocalStockBackend([db_paths["stock_synplanner"]])
    rxn_heldout = ReactionLookupBackend(db_paths["rxn_heldout_leakfree"])
    rxn_uspto_strict = ReactionLookupBackend(db_paths["rxn_uspto_strict_leakfree"])
    rxn_upper = ReactionLookupBackend(db_paths["rxn_retrieval_upper_bound"])
    rxn_uspto_single = ReactionLookupBackend(db_paths["rxn_uspto_single_step"])
    tb_17 = TemplateFallbackBackend(max_branch_per_node=16)

    lookup_fns = {
        "paroutes_n1": _make_offer_lookup(stock_n1),
        "paroutes_n5": _make_offer_lookup(stock_n5),
        "synplanner_sascore": _make_offer_lookup(stock_sp),
    }

    method_names = (
        "stock_only",
        "lookup_1step_smarts_2step",
        "retro_star_max4_leakage_free",
        "retro_star_max5_leakage_free",
        "retro_star_max5_uspto50k_only",
        "retro_star_max5_retrieval_upper_bound",
    )
    per_method_records: dict[str, list[dict[str, Any]]] = {m: [] for m in method_names}
    sc_weights = REPO_ROOT / DEFAULT_SCSCORE_WEIGHTS

    try:
        for t_info in all_targets:
            smi = t_info["smiles"]
            ik = t_info["inchikey"]
            ds = t_info["dataset"]
            offer_lookup = lookup_fns[ds]
            nm = normalize_smiles(smi)
            mol = nm.mol if nm is not None else None

            # 1. Stock-only lookup (max_steps=0)
            t0 = time.perf_counter()
            in_stock = _is_purchasable(smi, {}, offer_lookup, desalt_fallback=True)
            dt_stock = (time.perf_counter() - t0) * 1000.0
            per_method_records["stock_only"].append(
                {
                    **t_info,
                    "solved": in_stock,
                    "n_steps": 0 if in_stock else None,
                    "leaf_completeness": 1.0 if in_stock else 0.0,
                    "time_ms": dt_stock,
                }
            )

            # 2. Leakage-free 1-step lookup + 2-step greedy SMARTS (17 templates)
            t1 = time.perf_counter()
            if in_stock:
                per_method_records["lookup_1step_smarts_2step"].append(
                    {
                        **t_info,
                        "solved": True,
                        "n_steps": 0,
                        "leaf_completeness": 1.0,
                        "time_ms": (time.perf_counter() - t1) * 1000.0,
                    }
                )
            else:
                b_routes = _eval_baseline_1step_lookup_2step_smarts(
                    smi, ik, offer_lookup, rxn_heldout, tb_17
                )
                dt_b = (time.perf_counter() - t1) * 1000.0
                b_solved = [r for r in b_routes if r.leaves and all(r.leaves_purchasable)]
                if b_solved:
                    br = b_solved[0]
                    lcomp = sum(1 for x in br.leaves_purchasable if x) / len(br.leaves_purchasable)
                    per_method_records["lookup_1step_smarts_2step"].append(
                        {
                            **t_info,
                            "solved": True,
                            "n_steps": br.n_steps,
                            "leaf_completeness": lcomp,
                            "time_ms": dt_b,
                        }
                    )
                else:
                    best_partial = (
                        max(
                            sum(1 for x in r.leaves_purchasable if x) / len(r.leaves_purchasable)
                            for r in b_routes
                            if r.leaves_purchasable
                        )
                        if b_routes
                        else 0.0
                    )
                    per_method_records["lookup_1step_smarts_2step"].append(
                        {
                            **t_info,
                            "solved": False,
                            "n_steps": None,
                            "leaf_completeness": best_partial,
                            "time_ms": dt_b,
                        }
                    )

            # 3a, 3b (Headline Leakage-Free), 3c (Strict Zero-PaRoutes), 4 (Retrieval Upper Bound)
            for m_name, m_steps, m_rxn_db in (
                ("retro_star_max4_leakage_free", 4, rxn_heldout),
                ("retro_star_max5_leakage_free", 5, rxn_heldout),
                ("retro_star_max5_uspto50k_only", 5, rxn_uspto_strict),
                ("retro_star_max5_retrieval_upper_bound", 5, rxn_upper),
            ):
                ts = time.perf_counter()
                if in_stock:
                    per_method_records[m_name].append(
                        {
                            **t_info,
                            "solved": True,
                            "n_steps": 0,
                            "leaf_completeness": 1.0,
                            "time_ms": (time.perf_counter() - ts) * 1000.0,
                        }
                    )
                    continue
                rs_routes = find_routes(
                    smi,
                    ik,
                    max_steps=m_steps,
                    route_backend="auto",
                    planner="retro_star",
                    offer_lookup=offer_lookup,
                    lookup_backend=m_rxn_db,
                    template_backend=tb_17,
                    target_mol=mol,
                    scscore_weights=sc_weights,
                    retro_star_max_expansions=retro_star_expansions,
                )
                dt_rs = (time.perf_counter() - ts) * 1000.0
                rs_solved = [r for r in rs_routes if r.leaves and all(r.leaves_purchasable)]
                if rs_solved:
                    best_r = rs_solved[0]
                    lcomp = sum(1 for x in best_r.leaves_purchasable if x) / len(best_r.leaves_purchasable)
                    per_method_records[m_name].append(
                        {
                            **t_info,
                            "solved": True,
                            "n_steps": best_r.n_steps,
                            "leaf_completeness": lcomp,
                            "time_ms": dt_rs,
                        }
                    )
                else:
                    best_partial = (
                        max(
                            sum(1 for x in r.leaves_purchasable if x) / len(r.leaves_purchasable)
                            for r in rs_routes
                            if r.leaves_purchasable
                        )
                        if rs_routes
                        else 0.0
                    )
                    per_method_records[m_name].append(
                        {
                            **t_info,
                            "solved": False,
                            "n_steps": None,
                            "leaf_completeness": best_partial,
                            "time_ms": dt_rs,
                        }
                    )

        uspto50k_single_step = _evaluate_uspto50k_single_step_topk(
            corpora_dir / "syntheseus_uspto_50k" / "uspto50k.zip",
            rxn_uspto_single,
            n_test_samples=uspto50k_test_n,
            seed=seed,
        )
    finally:
        stock_n1.close()
        stock_n5.close()
        stock_sp.close()
        rxn_heldout.close()
        rxn_uspto_strict.close()
        rxn_upper.close()
        rxn_uspto_single.close()

    overall_summary = {m: _summarize_method_runs(per_method_records[m]) for m in method_names}

    by_dataset: dict[str, dict[str, Any]] = {}
    for ds_name in ("paroutes_n1", "paroutes_n5", "synplanner_sascore"):
        by_dataset[ds_name] = {
            m: _summarize_method_runs([r for r in per_method_records[m] if r["dataset"] == ds_name])
            for m in method_names
        }

    by_paroutes_depth: dict[str, dict[str, Any]] = {}
    for d in (2, 3, 4, 5):
        by_paroutes_depth[f"depth_{d}"] = {
            m: _summarize_method_runs([r for r in per_method_records[m] if r["ref_depth"] == d])
            for m in method_names
        }

    convergent_summary = {
        m: _summarize_method_runs([r for r in per_method_records[m] if r["is_convergent_ref"] is True])
        for m in method_names
    }

    by_sa_band: dict[str, dict[str, Any]] = {}
    for band in ("[1.5, 2.5)", "[2.5, 3.5)", "[3.5, 4.5)", "[4.5, 10.0]"):
        by_sa_band[band] = {
            m: _summarize_method_runs([r for r in per_method_records[m] if r["sa_band"] == band])
            for m in method_names
        }

    index_hygiene = {
        "stock_only": {
            "role": "baseline_stock_lookup",
            "reaction_index": "none",
            "target_reference_routes_in_index": False,
            "sampled_target_trees_excluded": excluded_trees_total,
            "heldout_paroutes_trees_indexed": 0,
            "indexed_reaction_rows": 0,
            "paroutes_target_product_inchikey_overlap": 0,
            "all_target_product_inchikey_overlap": 0,
            "all_target_connectivity_block_overlap": 0,
        },
        "lookup_1step_smarts_2step": {
            "role": "baseline_leakage_free",
            "reaction_index": "heldout_non_target_paroutes_plus_uspto50k_trainval",
            "target_reference_routes_in_index": False,
            "sampled_target_trees_excluded": excluded_trees_total,
            "heldout_paroutes_trees_indexed": heldout_trees_indexed,
            **hygiene_heldout,
            **filter_counts,
        },
        "retro_star_max4_leakage_free": {
            "role": "leakage_free_planner",
            "reaction_index": "heldout_non_target_paroutes_plus_uspto50k_trainval",
            "target_reference_routes_in_index": False,
            "sampled_target_trees_excluded": excluded_trees_total,
            "heldout_paroutes_trees_indexed": heldout_trees_indexed,
            **hygiene_heldout,
            **filter_counts,
        },
        "retro_star_max5_leakage_free": {
            "role": "headline_leakage_free_planner",
            "reaction_index": "heldout_non_target_paroutes_plus_uspto50k_trainval",
            "target_reference_routes_in_index": False,
            "sampled_target_trees_excluded": excluded_trees_total,
            "heldout_paroutes_trees_indexed": heldout_trees_indexed,
            **hygiene_heldout,
            **filter_counts,
        },
        "retro_star_max5_uspto50k_only": {
            "role": "strict_zero_paroutes_ablation",
            "reaction_index": "uspto50k_trainval_target_filtered_only",
            "target_reference_routes_in_index": False,
            "sampled_target_trees_excluded": excluded_trees_total,
            "heldout_paroutes_trees_indexed": 0,
            **hygiene_uspto,
        },
        "retro_star_max5_retrieval_upper_bound": {
            "role": "retrieval_upper_bound_not_headline",
            "reaction_index": "all_paroutes_trees_including_targets_plus_uspto50k_trainval",
            "target_reference_routes_in_index": True,
            "sampled_target_trees_excluded": 0,
            "heldout_paroutes_trees_indexed": total_trees_n1 + total_trees_n5,
            **hygiene_upper,
        },
    }

    return {
        "benchmark_schema": "buy-or-make/retrosynthesis-benchmark/v1",
        "seed": seed,
        "data_source": "open_stopping_stock",
        "environment": {
            "python_version": platform.python_version(),
            "rdkit_version": rdkit.__version__,
            "platform": platform.platform(),
            "processor": platform.processor() or platform.machine(),
        },
        "datasets": manifest,
        "multi_step_benchmark": {
            "headline_method": HEADLINE_METHOD,
            "n_core_smarts_templates": N_CORE_TEMPLATES,
            "n_extended_smarts_templates": N_EXT_TEMPLATES,
            "n_targets_total": len(all_targets),
            "target_counts_by_dataset": {
                "paroutes_n1": len(targets_n1),
                "paroutes_n5": len(targets_n5),
                "synplanner_sascore": len(targets_sp),
            },
            "index_hygiene": index_hygiene,
            "overall": overall_summary,
            "by_dataset": by_dataset,
            "by_paroutes_depth": by_paroutes_depth,
            "paroutes_convergent_subset": convergent_summary,
            "by_sascore_band": by_sa_band,
        },
        "uspto_50k_single_step": uspto50k_single_step,
    }


def render_markdown_benchmark_section(report: dict[str, Any]) -> str:
    """Render deterministic Markdown tables from `results/retrosynthesis_benchmark.json`."""
    ms = report["multi_step_benchmark"]
    n_tot = ms["n_targets_total"]
    n_ext = ms.get("n_extended_smarts_templates", N_EXT_TEMPLATES)
    ov = ms["overall"]
    by_ds = ms["by_dataset"]
    by_d = ms["by_paroutes_depth"]
    by_sa = ms["by_sascore_band"]
    hyg = ms["index_hygiene"]
    h_main = hyg["retro_star_max5_leakage_free"]
    u50 = report["uspto_50k_single_step"]
    u_m = u50["methods"]
    ds_meta = report["datasets"]

    method_labels = (
        ("stock_only", "1. Stock-only lookup (`max_steps=0`)"),
        (
            "lookup_1step_smarts_2step",
            "2. 1-step held-out lookup + 2-step greedy SMARTS (`17` templates, leakage-free)",
        ),
        (
            "retro_star_max4_leakage_free",
            f"3a. Leakage-free `Retro*` (`max_steps=4`, held-out `PaRoutes` + `USPTO-50K` + `{n_ext}` SMARTS)",
        ),
        (
            "retro_star_max5_leakage_free",
            f"3b. **Headline — Leakage-free `Retro*`** (`max_steps=5`, held-out `PaRoutes` + `USPTO-50K` + `{n_ext}` SMARTS)",
        ),
        (
            "retro_star_max5_uspto50k_only",
            f"3c. Strict zero-`PaRoutes` `Retro*` (`max_steps=5`, `USPTO-50K` train+val + `{n_ext}` SMARTS only)",
        ),
        (
            "retro_star_max5_retrieval_upper_bound",
            "4. *Reference-route retrieval upper bound* (`max_steps=5`, target's own `PaRoutes` route indexed — NOT headline)",
        ),
    )

    lines: list[str] = [
        BENCHMARK_START_MARKER,
        (
            f"Measured offline on CPU (`Python {report['environment']['python_version']}`, "
            f"`RDKit {report['environment']['rdkit_version']}`, seed `{report['seed']}`) via "
            f"[`scripts/benchmark_retrosynthesis.py`](scripts/benchmark_retrosynthesis.py) -> "
            f"[`results/retrosynthesis_benchmark.json`](results/retrosynthesis_benchmark.json) over "
            f"**{n_tot} stratified multi-step targets** (`PaRoutes 2.0 n1` `n={ms['target_counts_by_dataset']['paroutes_n1']}`, "
            f"`PaRoutes 2.0 n5` `n={ms['target_counts_by_dataset']['paroutes_n5']}`, "
            f"`SynPlanner` SAscore-stratified `n={ms['target_counts_by_dataset']['synplanner_sascore']}`) plus "
            f"**{u50['n_test_targets']} held-out `Syntheseus USPTO-50K` test reactions** (`raw_test.csv` strictly excluded). "
            f"**Split & target-exclusion hygiene**: the headline leakage-free reaction index (`{h_main['indexed_reaction_rows']:,}` reactions "
            f"across `{h_main['heldout_paroutes_trees_indexed']:,}` non-target `PaRoutes` trees + `USPTO-50K` train+val) excludes all "
            f"`{h_main['sampled_target_trees_excluded']}` sampled `PaRoutes` target trees **and** any reaction whose product `InChIKey` "
            f"or 14-char connectivity block matches any of the `{n_tot}` evaluated targets "
            f"(`all_target_product_inchikey_overlap = {h_main['all_target_product_inchikey_overlap']}`):"
        ),
        "",
        "#### Overall & Dataset-Stratified Multi-Step Solve Rates (`n = " + str(n_tot) + "`)",
        "",
        "| Planner / Configuration | Target Routes Indexed? | Overall Solve % | `PaRoutes n1` | `PaRoutes n5` | `SynPlanner` | Leaf Stock Completeness | Median Search Time |",
        "|---|:---:|---:|---:|---:|---:|---:|---:|",
    ]

    for key, label in method_labels:
        m_ov = ov[key]
        m_n1 = by_ds["paroutes_n1"][key]
        m_n5 = by_ds["paroutes_n5"][key]
        m_sp = by_ds["synplanner_sascore"][key]
        in_idx = "Yes (upper bound)" if hyg[key]["target_reference_routes_in_index"] else "**No (0 overlap)**"
        lines.append(
            f"| {label} | {in_idx} | {100.0 * m_ov['solve_rate']:.1f}% ({m_ov['solved_count']}/{m_ov['n_targets']}) | "
            f"{100.0 * m_n1['solve_rate']:.1f}% | {100.0 * m_n5['solve_rate']:.1f}% | "
            f"{100.0 * m_sp['solve_rate']:.1f}% | "
            f"{100.0 * m_ov['mean_leaf_stock_completeness']:.1f}% (solved: {100.0 * m_ov['solved_leaf_stock_completeness']:.1f}%) | "
            f"{m_ov['median_time_ms']:.2f} ms |"
        )

    lines.extend(
        [
            "",
            "#### Solve Rate Stratified by `PaRoutes 2.0` Reference Route Depth & `SAscore` Band",
            "",
            "| Planner | Depth 2 | Depth 3 | Depth 4 | Depth 5 | Convergent Subset | `SA [1.5,2.5)` | `SA [2.5,3.5)` | `SA [3.5,4.5)` | `SA [4.5,10]` |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    conv = ms["paroutes_convergent_subset"]
    for key, label in method_labels:
        lines.append(
            f"| {label} | "
            f"{100.0 * by_d['depth_2'][key]['solve_rate']:.1f}% | "
            f"{100.0 * by_d['depth_3'][key]['solve_rate']:.1f}% | "
            f"{100.0 * by_d['depth_4'][key]['solve_rate']:.1f}% | "
            f"{100.0 * by_d['depth_5'][key]['solve_rate']:.1f}% | "
            f"{100.0 * conv[key]['solve_rate']:.1f}% | "
            f"{100.0 * by_sa['[1.5, 2.5)'][key]['solve_rate']:.1f}% | "
            f"{100.0 * by_sa['[2.5, 3.5)'][key]['solve_rate']:.1f}% | "
            f"{100.0 * by_sa['[3.5, 4.5)'][key]['solve_rate']:.1f}% | "
            f"{100.0 * by_sa['[4.5, 10.0]'][key]['solve_rate']:.1f}% |"
        )

    ext_u_key = f"retro_star_lookup_trainval_plus_{u50.get('n_extended_smarts_templates', n_ext)}_smarts"
    b_u = u_m["lookup_trainval_plus_17_smarts"]
    e_u = u_m[ext_u_key]
    u_n_ext = u50.get("n_extended_smarts_templates", n_ext)
    lines.extend(
        [
            "",
            f"#### Held-Out `Syntheseus USPTO-50K` Single-Step Exact Reactant-Set Accuracy (`n = {u50['n_test_targets']}`)",
            "",
            "| Single-Step Disconnection Engine | Top-1 | Top-3 | Top-5 | Top-10 | Median Time |",
            "|---|---:|---:|---:|---:|---:|",
            f"| `USPTO-50K` train+val lookup + `17` core SMARTS | {100.0 * b_u['top_1']:.1f}% | {100.0 * b_u['top_3']:.1f}% | {100.0 * b_u['top_5']:.1f}% | {100.0 * b_u['top_10']:.1f}% | {b_u['median_time_ms']:.2f} ms |",
            f"| `USPTO-50K` train+val lookup + `{u_n_ext}` `Retro*` SMARTS | {100.0 * e_u['top_1']:.1f}% | {100.0 * e_u['top_3']:.1f}% | {100.0 * e_u['top_5']:.1f}% | {100.0 * e_u['top_10']:.1f}% | {e_u['median_time_ms']:.2f} ms |",
            "",
            "#### Verified Input Dataset SHA-256 Digests",
            "",
            "| Artifact | License | Bytes | SHA-256 Digest |",
            "|---|---|---:|---|",
        ]
    )
    for k in (
        "paroutes_stock_n1",
        "paroutes_stock_n5",
        "paroutes_routes_n1",
        "paroutes_routes_n5",
        "syntheseus_uspto_50k",
        "synplanner_building_blocks",
        "synplanner_targets_sa_1_5_2_5",
        "synplanner_targets_sa_2_5_3_5",
        "synplanner_targets_sa_3_5_4_5",
        "synplanner_targets_sa_4_5_5_5",
        "scscore_weights",
    ):
        info = ds_meta[k]
        lines.append(
            f"| `{info['filename']}` (`{k}`) | `{info['spdx_license']}` | {info['size_bytes']:,} | `{info['sha256']}` |"
        )
    lines.append(BENCHMARK_END_MARKER)
    return "\n".join(lines)


def update_readme_benchmark_block(
    report: dict[str, Any],
    readme_path: Path = README_PATH,
) -> None:
    rendered = render_markdown_benchmark_section(report)
    text = readme_path.read_text(encoding="utf-8")
    if BENCHMARK_START_MARKER in text and BENCHMARK_END_MARKER in text:
        before = text.split(BENCHMARK_START_MARKER, 1)[0]
        after = text.split(BENCHMARK_END_MARKER, 1)[1]
        updated = before + rendered + after
    else:
        anchor = "## Demo coverage and CPU runtime"
        section = (
            "## Multi-Step Retrosynthesis & Single-Step Benchmark (`PaRoutes 2.0`, `SynPlanner`, `USPTO-50K`)\n\n"
            + rendered
            + "\n\n"
        )
        if anchor in text:
            updated = text.replace(anchor, section + anchor, 1)
        else:
            updated = text.rstrip() + "\n\n" + section
    readme_path.write_text(updated, encoding="utf-8")


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run reproducible leakage-free multi-step & single-step retrosynthesis benchmark."
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_RESULTS_JSON,
        help="Destination JSON path (default: results/retrosynthesis_benchmark.json).",
    )
    parser.add_argument(
        "--corpora-dir",
        type=Path,
        default=DEFAULT_OPEN_CORPORA_DIR,
        help="Directory for gitignored raw corpora (default: data/open_corpora).",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("data/benchmark_cache"),
        help="Directory for gitignored SQLite benchmark indexes.",
    )
    parser.add_argument(
        "--update-readme",
        action="store_true",
        help="Update README.md benchmark tables from the generated JSON.",
    )
    args = parser.parse_args(argv)

    report = run_full_benchmark(corpora_dir=args.corpora_dir, cache_dir=args.cache_dir)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote benchmark report -> {args.out}")
    if args.update_readme:
        update_readme_benchmark_block(report, README_PATH)
        print(f"Updated benchmark tables in {README_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
