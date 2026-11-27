#!/usr/bin/env python3
"""Reproducible CPU multi-step & single-step retrosynthesis benchmark on real public corpora.

Evaluates `buy-or-make` across:
1. `PaRoutes 2.0` (`n1` and `n5` stopping stocks + multi-step patent route targets, Zenodo `7341155`, `CC-BY-4.0`)
2. `SynPlanner` SAscore-stratified retrosynthesis benchmark targets + building blocks (`MIT`)
3. `Syntheseus USPTO-50K` held-out `raw_test.csv` single-step top-k accuracy (`MIT`, train+val only indexed)

Compares:
- `stock_only`: Direct stopping-stock lookup (`max_steps=0`)
- `lookup_1step_smarts_2step`: 1-step reaction lookup + 2-step greedy SMARTS (`17` templates)
- `retro_star_max4`: Chen et al. `Retro*` best-first AND-OR search (`max_steps=4`) over indexed
  `PaRoutes`/`USPTO-50K` + `44` atom-conserved SMARTS templates
- `retro_star_max5`: Chen et al. `Retro*` best-first AND-OR search (`max_steps=5`) over indexed
  `PaRoutes`/`USPTO-50K` + `44` atom-conserved SMARTS templates
- `retro_star_max5_uspto50k_only`: Zero-target-route ablation (`Retro*`, `max_steps=5` using ONLY
  `USPTO-50K` train+val + `44` SMARTS templates, holding out all `PaRoutes` reference reactions)

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
from typing import Any, Callable, Optional

import rdkit

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from buyormake.canonical import file_sha256, normalize_smiles  # noqa: E402
from buyormake.open_corpora import (  # noqa: E402
    CORPUS_ARTIFACTS,
    DEFAULT_OPEN_CORPORA_DIR,
    _index_reaction_pairs,
    download_corpus_artifact,
    import_paroutes_stock,
    import_synplanner_building_blocks,
    load_synplanner_targets,
    parse_paroutes_route_tree,
)
from buyormake.retro_star import RetroStarPlanner  # noqa: E402
from buyormake.routes import (  # noqa: E402
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


def _sample_paroutes_targets(
    routes_json_path: Path,
    *,
    dataset_label: str,
    per_depth: int,
    depths: tuple[int, ...] = (2, 3, 4, 5),
    seed: int = 20261127,
) -> tuple[list[dict[str, Any]], list[tuple[str, str]]]:
    """Deterministically sample stratified multi-step patent targets by route depth `{2,3,4,5}`
    and extract deduplicated `(rxn_smiles, ref)` pairs across the PaRoutes corpus.
    """
    raw_trees = json.loads(routes_json_path.read_text(encoding="utf-8"))
    by_depth: dict[int, list[dict[str, Any]]] = {d: [] for d in depths}
    seen_target_iks: set[str] = set()
    rxn_pairs_seen: set[tuple[str, str]] = set()
    rxn_pairs: list[tuple[str, str]] = []

    def _collect_tree_reactions(node: dict[str, Any], prefix: str) -> None:
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
                        if (rsmi, ref) not in rxn_pairs_seen:
                            rxn_pairs_seen.add((rsmi, ref))
                            rxn_pairs.append((rsmi, ref))
                    for mc in mol_ch:
                        _collect_tree_reactions(mc, prefix)

    for idx, item in enumerate(raw_trees):
        root = item.get("dict") if isinstance(item.get("dict"), dict) else item
        if isinstance(root, dict):
            _collect_tree_reactions(root, f"{dataset_label}:r{idx}")

    # Deterministically shuffle indices to avoid file-ordering bias before parsing stratified targets
    rng = random.Random(seed + (101 if "n1" in dataset_label else 202))
    indices = list(range(len(raw_trees)))
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
        by_depth[route.n_steps].append(
            {
                "smiles": nm.canonical_smiles,
                "inchikey": nm.inchikey,
                "dataset": dataset_label,
                "ref_depth": route.n_steps,
                "sascore": round(float(sa_val), 4),
                "sa_band": _sa_band_label(sa_val),
                "is_convergent_ref": _is_convergent_route(route),
            }
        )

    targets: list[dict[str, Any]] = []
    for d in depths:
        targets.extend(by_depth[d])
    return targets, rxn_pairs


def _sa_band_label(sa: float) -> str:
    if sa < 2.5:
        return "[1.5, 2.5)"
    if sa < 3.5:
        return "[2.5, 3.5)"
    if sa < 4.5:
        return "[3.5, 4.5)"
    return "[4.5, 10.0]"


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


def _build_benchmark_databases(
    corpora_dir: Path,
    cache_dir: Path,
    paroutes_rxn_pairs: list[tuple[str, str]],
) -> dict[str, Path]:
    """Build gitignored SQLite stopping-stock and split-clean reaction indexes."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    stock_n1_db = cache_dir / "bench_stock_n1.sqlite"
    stock_n5_db = cache_dir / "bench_stock_n5.sqlite"
    stock_sp_db = cache_dir / "bench_stock_synplanner.sqlite"
    rxn_uspto_trainval_db = cache_dir / "bench_rxn_uspto50k_trainval.sqlite"
    rxn_combined_db = cache_dir / "bench_rxn_paroutes_uspto50k.sqlite"

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

    uspto_zip = corpora_dir / "syntheseus_uspto_50k" / "uspto50k.zip"
    if not rxn_uspto_trainval_db.exists():
        def _iter_trainval():
            with zipfile.ZipFile(uspto_zip, "r") as zf:
                for split_csv in ("raw_train.csv", "raw_val.csv"):
                    with zf.open(split_csv, "r") as raw:
                        rdr = csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"))
                        for i, row in enumerate(rdr):
                            rxn = row.get("reactants>reagents>production")
                            ref = row.get("id") or f"{split_csv}:{i}"
                            yield rxn, ref

        _index_reaction_pairs(
            _iter_trainval(),
            db_path=rxn_uspto_trainval_db,
            step_source="syntheseus_uspto_50k_lookup",
            default_prefix="uspto50k",
        )

    if not rxn_combined_db.exists():
        # Copy train+val USPTO-50K reactions and add indexed PaRoutes patent reactions
        conn_src = sqlite3.connect(str(rxn_uspto_trainval_db))
        conn_dst = sqlite3.connect(str(rxn_combined_db))
        conn_src.backup(conn_dst)
        conn_dst.close()
        conn_src.close()
        _index_reaction_pairs(
            paroutes_rxn_pairs,
            db_path=rxn_combined_db,
            step_source="paroutes_lookup",
            default_prefix="paroutes",
        )

    return {
        "stock_n1": stock_n1_db,
        "stock_n5": stock_n5_db,
        "stock_synplanner": stock_sp_db,
        "rxn_uspto_trainval": rxn_uspto_trainval_db,
        "rxn_combined": rxn_combined_db,
    }


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
    """Method 2: 1-step reaction lookup + 2-step greedy SMARTS (`17` templates)."""
    offer_cache: dict[str, list[Offer]] = {}
    routes: list[Route] = []

    # 1. strictly 1-step reaction lookup
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

    # 2. 2-step greedy SMARTS fallback (`17` templates)
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
                rxn = row.get("reactants>reagents>production")
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

    tb_17 = TemplateFallbackBackend(max_branch_per_node=25)
    rs_ext = RetroStarPlanner(max_branch_per_node=25)

    ks = (1, 3, 5, 10)
    base_hits = {k: 0 for k in ks}
    ext_hits = {k: 0 for k in ks}
    base_times: list[float] = []
    ext_times: list[float] = []

    for smi, ik, mol, gt_set in test_items:
        # 1. Baseline: train+val lookup + 17 RETRO_TEMPLATES
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

        # 2. Extended: train+val lookup + 44 Retro* atom-conserved medicinal chemistry SMARTS
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
    return {
        "dataset": "syntheseus_uspto_50k_test",
        "split_hygiene": "raw_train.csv + raw_val.csv indexed; raw_test.csv strictly held out (0 test leakage)",
        "n_test_targets": total,
        "methods": {
            "lookup_trainval_plus_17_smarts": {
                "top_1": round(base_hits[1] / total, 4),
                "top_3": round(base_hits[3] / total, 4),
                "top_5": round(base_hits[5] / total, 4),
                "top_10": round(base_hits[10] / total, 4),
                "median_time_ms": round(float(statistics.median(base_times)), 3),
            },
            "retro_star_lookup_trainval_plus_44_smarts": {
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
    seed: int = 20261127,
) -> dict[str, Any]:
    manifest = _ensure_artifacts(corpora_dir)

    targets_n1, rxn_n1 = _sample_paroutes_targets(
        corpora_dir / "paroutes" / "ref_routes_n1.json",
        dataset_label="paroutes_n1",
        per_depth=paroutes_per_depth,
        seed=seed,
    )
    targets_n5, rxn_n5 = _sample_paroutes_targets(
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

    db_paths = _build_benchmark_databases(
        corpora_dir, cache_dir, [*rxn_n1, *rxn_n5]
    )

    stock_n1 = LocalStockBackend([db_paths["stock_n1"]])
    stock_n5 = LocalStockBackend([db_paths["stock_n5"]])
    stock_sp = LocalStockBackend([db_paths["stock_synplanner"]])
    rxn_combined = ReactionLookupBackend(db_paths["rxn_combined"])
    rxn_trainval = ReactionLookupBackend(db_paths["rxn_uspto_trainval"])
    tb_17 = TemplateFallbackBackend(max_branch_per_node=16)

    lookup_fns = {
        "paroutes_n1": _make_offer_lookup(stock_n1),
        "paroutes_n5": _make_offer_lookup(stock_n5),
        "synplanner_sascore": _make_offer_lookup(stock_sp),
    }

    method_names = (
        "stock_only",
        "lookup_1step_smarts_2step",
        "retro_star_max4",
        "retro_star_max5",
        "retro_star_max5_uspto50k_only",
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

            # 2. 1-step reaction lookup + 2-step greedy SMARTS (17 templates)
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
                    smi, ik, offer_lookup, rxn_combined, tb_17
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

            # 3, 4, 5. Retro* (max_steps=4, max_steps=5 over PaRoutes+USPTO50K+SMARTS, and max_steps=5 USPTO50K+SMARTS only)
            for m_name, m_steps, m_rxn_db in (
                ("retro_star_max4", 4, rxn_combined),
                ("retro_star_max5", 5, rxn_combined),
                ("retro_star_max5_uspto50k_only", 5, rxn_trainval),
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
                    retro_star_max_expansions=80,
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
            rxn_trainval,
            n_test_samples=uspto50k_test_n,
            seed=seed,
        )
    finally:
        stock_n1.close()
        stock_n5.close()
        stock_sp.close()
        rxn_combined.close()
        rxn_trainval.close()

    # Build stratified summaries
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
            "n_targets_total": len(all_targets),
            "target_counts_by_dataset": {
                "paroutes_n1": len(targets_n1),
                "paroutes_n5": len(targets_n5),
                "synplanner_sascore": len(targets_sp),
            },
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
    ov = ms["overall"]
    by_ds = ms["by_dataset"]
    by_d = ms["by_paroutes_depth"]
    by_sa = ms["by_sascore_band"]
    u50 = report["uspto_50k_single_step"]
    u_m = u50["methods"]
    ds_meta = report["datasets"]

    method_labels = (
        ("stock_only", "1. Stock-only lookup (`max_steps=0`)"),
        ("lookup_1step_smarts_2step", "2. 1-step lookup + 2-step greedy SMARTS (`17` templates)"),
        ("retro_star_max4", "3a. `Retro*` AND-OR search (`max_steps=4`, `PaRoutes`+`USPTO-50K`+`44` SMARTS)"),
        ("retro_star_max5", "3b. `Retro*` AND-OR search (`max_steps=5`, `PaRoutes`+`USPTO-50K`+`44` SMARTS)"),
        ("retro_star_max5_uspto50k_only", "3c. `Retro*` zero-target-route ablation (`max_steps=5`, `USPTO-50K`+`44` SMARTS only)"),
    )

    lines: list[str] = [
        BENCHMARK_START_MARKER,
        f"Measured offline on CPU (`Python {report['environment']['python_version']}`, `RDKit {report['environment']['rdkit_version']}`, seed `{report['seed']}`) via [`scripts/benchmark_retrosynthesis.py`](scripts/benchmark_retrosynthesis.py) -> [`results/retrosynthesis_benchmark.json`](results/retrosynthesis_benchmark.json) over **{n_tot} stratified multi-step targets** (`PaRoutes 2.0 n1` `n={ms['target_counts_by_dataset']['paroutes_n1']}`, `PaRoutes 2.0 n5` `n={ms['target_counts_by_dataset']['paroutes_n5']}`, `SynPlanner` SAscore-stratified `n={ms['target_counts_by_dataset']['synplanner_sascore']}`) plus **{u50['n_test_targets']} held-out `Syntheseus USPTO-50K` test reactions** (`raw_test.csv` strictly excluded from SQLite reaction indexing):",
        "",
        "#### Overall & Dataset-Stratified Multi-Step Solve Rates (`n = " + str(n_tot) + "`)",
        "",
        "| Planner / Configuration | Overall Solve % | `PaRoutes n1` | `PaRoutes n5` | `SynPlanner` | Leaf Stock Completeness | Median Search Time |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]

    for key, label in method_labels:
        m_ov = ov[key]
        m_n1 = by_ds["paroutes_n1"][key]
        m_n5 = by_ds["paroutes_n5"][key]
        m_sp = by_ds["synplanner_sascore"][key]
        lines.append(
            f"| {label} | {100.0 * m_ov['solve_rate']:.1f}% ({m_ov['solved_count']}/{m_ov['n_targets']}) | "
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

    b_u = u_m["lookup_trainval_plus_17_smarts"]
    e_u = u_m["retro_star_lookup_trainval_plus_44_smarts"]
    lines.extend(
        [
            "",
            f"#### Held-Out `Syntheseus USPTO-50K` Single-Step Exact Reactant-Set Accuracy (`n = {u50['n_test_targets']}`)",
            "",
            "| Single-Step Disconnection Engine | Top-1 | Top-3 | Top-5 | Top-10 | Median Time |",
            "|---|---:|---:|---:|---:|---:|",
            f"| `USPTO-50K` train+val lookup + `17` core SMARTS | {100.0 * b_u['top_1']:.1f}% | {100.0 * b_u['top_3']:.1f}% | {100.0 * b_u['top_5']:.1f}% | {100.0 * b_u['top_10']:.1f}% | {b_u['median_time_ms']:.2f} ms |",
            f"| `USPTO-50K` train+val lookup + `44` `Retro*` SMARTS | {100.0 * e_u['top_1']:.1f}% | {100.0 * e_u['top_3']:.1f}% | {100.0 * e_u['top_5']:.1f}% | {100.0 * e_u['top_10']:.1f}% | {e_u['median_time_ms']:.2f} ms |",
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
        description="Run reproducible multi-step & single-step retrosynthesis benchmark."
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
