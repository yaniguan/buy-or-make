"""Regression guard verifying zero drift between `results/retrosynthesis_benchmark.json` and `README.md`,
strict multi-step reaction-index leakage exclusion (Issue #35), all 19/19 pinned corpus digests,
and repository-wide pure-English (zero CJK) hygiene checks.
"""

from __future__ import annotations

import csv
import importlib.util
import io
import json
import re
import zipfile
from pathlib import Path

from buyormake.canonical import normalize_smiles
from buyormake.open_corpora import CORPUS_ARTIFACTS
from buyormake.sources import SOURCES

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_JSON = REPO_ROOT / "results" / "retrosynthesis_benchmark.json"
README_MD = REPO_ROOT / "README.md"
BENCHMARK_SCRIPT = REPO_ROOT / "scripts" / "benchmark_retrosynthesis.py"


def _load_benchmark_script():
    spec = importlib.util.spec_from_file_location("benchmark_retrosynthesis", BENCHMARK_SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _has_cjk(text: str) -> bool:
    return any(
        "\u4e00" <= ch <= "\u9fff"
        or "\u3400" <= ch <= "\u4dbf"
        or "\u3040" <= ch <= "\u30ff"
        or "\uac00" <= ch <= "\ud7af"
        for ch in text
    )


def test_all_19_corpus_artifacts_have_pinned_cryptographic_digest() -> None:
    """Issue #35: all 19/19 registered corpus artifacts carry a pinned expected_sha256 or expected_md5."""
    assert len(CORPUS_ARTIFACTS) == 19
    sha_re = re.compile(r"^[0-9a-f]{64}$")
    md5_re = re.compile(r"^[0-9a-f]{32}$")
    for key, art in CORPUS_ARTIFACTS.items():
        has_sha = art.expected_sha256 is not None and bool(sha_re.match(art.expected_sha256))
        has_md5 = art.expected_md5 is not None and bool(md5_re.match(art.expected_md5))
        assert has_sha or has_md5, (
            f"Corpus artifact {key!r} lacks a valid pinned expected_sha256 or expected_md5"
        )


def test_headline_reaction_index_excludes_evaluated_target_trees_and_target_inchikeys(
    tmp_path: Path,
) -> None:
    """Issue #35: verify zero overlap between headline indexed reaction product InChIKeys and evaluated targets."""
    mod = _load_benchmark_script()

    corpora_dir = tmp_path / "open_corpora"
    (corpora_dir / "paroutes").mkdir(parents=True)
    (corpora_dir / "synplanner_data").mkdir(parents=True)
    (corpora_dir / "syntheseus_uspto_50k").mkdir(parents=True)

    # Target T1 (2-step route in tree 0): N-(4-(pyridin-2-yl)phenyl)acetamide
    t1_smi = "CC(=O)Nc1ccc(-c2ccccn2)cc1"
    i1_smi = "Nc1ccc(-c2ccccn2)cc1"
    # Non-target held-out tree 1: contains both a genuine non-target reaction AND a duplicate reaction producing T1!
    non_target_prod = "O=C(Nc1ccc(F)cc1)c1ccccc1"

    tree_target = {
        "type": "mol",
        "smiles": t1_smi,
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "metadata": {"ID": "rxn_t1_step1"},
                "children": [
                    {"type": "mol", "smiles": "CC(=O)O", "in_stock": True},
                    {
                        "type": "mol",
                        "smiles": i1_smi,
                        "in_stock": False,
                        "children": [
                            {
                                "type": "reaction",
                                "metadata": {"ID": "rxn_t1_step2"},
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
    tree_heldout = {
        "type": "mol",
        "smiles": non_target_prod,
        "in_stock": False,
        "children": [
            {
                "type": "reaction",
                "metadata": {"ID": "rxn_heldout_ok"},
                "children": [
                    {"type": "mol", "smiles": "O=C(O)c1ccccc1", "in_stock": True},
                    {
                        "type": "mol",
                        "smiles": t1_smi,
                        "in_stock": False,
                        "children": [
                            {
                                "type": "reaction",
                                "metadata": {"ID": "rxn_cross_tree_leak_trap"},
                                "children": [
                                    {"type": "mol", "smiles": "CC(=O)Cl", "in_stock": True},
                                    {"type": "mol", "smiles": i1_smi, "in_stock": True},
                                ],
                            }
                        ],
                    },
                ],
            }
        ],
    }

    for sub in ("n1", "n5"):
        (corpora_dir / "paroutes" / f"stock_{sub}.txt").write_text(
            "CC(=O)O\nNc1ccc(Br)cc1\nOB(O)c1ccccn1\nO=C(O)c1ccccc1\nNc1ccc(F)cc1\n",
            encoding="utf-8",
        )
        (corpora_dir / "paroutes" / f"ref_routes_{sub}.json").write_text(
            json.dumps([tree_target, tree_heldout]),
            encoding="utf-8",
        )

    (corpora_dir / "synplanner_data" / "building_blocks_em_sa_ln.smi").write_text(
        "CC(=O)O\nNc1ccccc1\n", encoding="utf-8"
    )

    # Also inject a USPTO-50K train reaction whose product equals T1 (cross-corpus leak trap!)
    uspto_zip = corpora_dir / "syntheseus_uspto_50k" / "uspto50k.zip"
    with zipfile.ZipFile(uspto_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for split_name in ("raw_train.csv", "raw_val.csv"):
            buf = io.StringIO()
            writer = csv.writer(buf)
            writer.writerow(["id", "reactants>reagents>production"])
            writer.writerow(["u_leak", f"CC(=O)O.{i1_smi}>>{t1_smi}"])
            writer.writerow(["u_clean", "O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1"])
            zf.writestr(split_name, buf.getvalue())

    targets_n1, tagged_n1, sampled_n1, _ = mod._sample_paroutes_targets(
        corpora_dir / "paroutes" / "ref_routes_n1.json",
        dataset_label="paroutes_n1",
        per_depth=1,
        depths=(2,),
        seed=20261127,
    )
    assert len(targets_n1) == 1
    t1_ik = normalize_smiles(t1_smi).inchikey  # type: ignore[union-attr]
    assert targets_n1[0]["inchikey"] == t1_ik
    assert sampled_n1 == {0}

    db_paths, _ = mod._build_benchmark_databases(
        corpora_dir,
        tmp_path / "cache",
        tagged_rxn_n1=tagged_n1,
        sampled_trees_n1=sampled_n1,
        tagged_rxn_n5=[],
        sampled_trees_n5=set(),
        evaluated_target_inchikeys={t1_ik},
    )

    hyg_headline = mod.verify_reaction_index_hygiene(
        db_paths["rxn_heldout_leakfree"],
        paroutes_target_inchikeys={t1_ik},
        all_target_inchikeys={t1_ik},
    )
    assert hyg_headline["indexed_reaction_rows"] >= 2
    assert hyg_headline["paroutes_target_product_inchikey_overlap"] == 0
    assert hyg_headline["all_target_product_inchikey_overlap"] == 0
    assert hyg_headline["all_target_connectivity_block_overlap"] == 0


def test_retrosynthesis_benchmark_json_schema_and_metrics() -> None:
    assert BENCHMARK_JSON.exists(), (
        "Missing committed results/retrosynthesis_benchmark.json; run `python scripts/benchmark_retrosynthesis.py --update-readme`"
    )
    report = json.loads(BENCHMARK_JSON.read_text(encoding="utf-8"))
    assert report["benchmark_schema"] == "buy-or-make/retrosynthesis-benchmark/v1"
    assert report["data_source"] == "open_stopping_stock"

    datasets = report["datasets"]
    sha_re = re.compile(r"^[0-9a-f]{64}$")
    for required_key in (
        "paroutes_stock_n1",
        "paroutes_stock_n5",
        "paroutes_routes_n1",
        "paroutes_routes_n5",
        "syntheseus_uspto_50k",
        "synplanner_building_blocks",
        "synplanner_targets_sa_1_5_2_5",
        "scscore_weights",
    ):
        assert required_key in datasets
        entry = datasets[required_key]
        assert sha_re.match(entry["sha256"]), f"Invalid SHA-256 on {required_key}: {entry['sha256']}"
        assert entry["size_bytes"] > 0

    ms = report["multi_step_benchmark"]
    assert ms["n_targets_total"] >= 200
    assert ms["headline_method"] == "retro_star_max5_leakage_free"

    # Verify index_hygiene is recorded for every method and headline has 0 target overlap
    hygiene = ms["index_hygiene"]
    headline_hyg = hygiene["retro_star_max5_leakage_free"]
    assert headline_hyg["target_reference_routes_in_index"] is False
    assert headline_hyg["paroutes_target_product_inchikey_overlap"] == 0
    assert headline_hyg["all_target_product_inchikey_overlap"] == 0
    assert headline_hyg["all_target_connectivity_block_overlap"] == 0
    assert headline_hyg["sampled_target_trees_excluded"] == 160
    assert hygiene["retro_star_max5_retrieval_upper_bound"]["target_reference_routes_in_index"] is True

    overall = ms["overall"]
    for m_key in (
        "stock_only",
        "lookup_1step_smarts_2step",
        "retro_star_max4_leakage_free",
        "retro_star_max5_leakage_free",
        "retro_star_max5_uspto50k_only",
        "retro_star_max5_retrieval_upper_bound",
    ):
        assert m_key in overall
        assert m_key in hygiene
        assert overall[m_key]["n_targets"] == ms["n_targets_total"]
        assert 0.0 <= overall[m_key]["solve_rate"] <= 1.0
        assert overall[m_key]["median_time_ms"] >= 0.0
        assert overall[m_key]["solved_leaf_stock_completeness"] == 1.0

    # Leakage-free headline monotonicity assertions against baselines and upper bound
    assert (
        overall["retro_star_max5_retrieval_upper_bound"]["solve_rate"]
        >= overall["retro_star_max5_leakage_free"]["solve_rate"]
    )
    assert (
        overall["retro_star_max5_leakage_free"]["solve_rate"]
        >= overall["retro_star_max4_leakage_free"]["solve_rate"]
    )
    assert (
        overall["retro_star_max5_leakage_free"]["solve_rate"]
        >= overall["retro_star_max5_uspto50k_only"]["solve_rate"]
    )
    assert (
        overall["retro_star_max5_leakage_free"]["solve_rate"]
        > overall["lookup_1step_smarts_2step"]["solve_rate"]
    )
    assert overall["lookup_1step_smarts_2step"]["solve_rate"] > overall["stock_only"]["solve_rate"]

    # Held-out USPTO-50K single-step top-k accuracy monotonicity
    u50 = report["uspto_50k_single_step"]
    assert u50["n_test_targets"] >= 200
    n_ext = u50["n_extended_smarts_templates"]
    ext_m = u50["methods"][f"retro_star_lookup_trainval_plus_{n_ext}_smarts"]
    base_m = u50["methods"]["lookup_trainval_plus_17_smarts"]
    assert ext_m["top_1"] <= ext_m["top_3"] <= ext_m["top_5"] <= ext_m["top_10"]
    assert ext_m["top_10"] >= base_m["top_10"]


def test_readme_benchmark_table_has_zero_drift_against_json() -> None:
    report = json.loads(BENCHMARK_JSON.read_text(encoding="utf-8"))
    mod = _load_benchmark_script()
    expected_block = mod.render_markdown_benchmark_section(report)

    readme_text = README_MD.read_text(encoding="utf-8")
    assert mod.BENCHMARK_START_MARKER in readme_text, "Missing BENCHMARK_TABLE_START in README.md"
    assert mod.BENCHMARK_END_MARKER in readme_text, "Missing BENCHMARK_TABLE_END in README.md"

    start_idx = readme_text.index(mod.BENCHMARK_START_MARKER)
    end_idx = readme_text.index(mod.BENCHMARK_END_MARKER) + len(mod.BENCHMARK_END_MARKER)
    actual_block = readme_text[start_idx:end_idx]

    assert actual_block == expected_block, (
        "README.md benchmark tables drifted from results/retrosynthesis_benchmark.json; "
        "re-run `python scripts/benchmark_retrosynthesis.py --update-readme`"
    )


def test_zero_cjk_characters_and_source_count_sync() -> None:
    readme_text = README_MD.read_text(encoding="utf-8")
    assert "| PubChem PUG-REST & PUG-View |" in readme_text
    assert len(SOURCES) == 19

    scan_dirs = (
        "buyormake",
        "community_extensions",
        "docs",
        "examples",
        "lean",
        "results",
        "schema",
        "scripts",
        "tests",
    )
    checked_files = [
        REPO_ROOT / "README.md",
        REPO_ROOT / "autolab_wire.py",
        REPO_ROOT / "pyproject.toml",
        REPO_ROOT / ".gitignore",
    ]
    for dname in scan_dirs:
        dpath = REPO_ROOT / dname
        if not dpath.exists():
            continue
        for p in dpath.rglob("*"):
            if p.is_file() and ".lake" not in p.parts and "__pycache__" not in p.parts:
                if p.suffix in {".py", ".md", ".json", ".lean", ".toml", ".yml", ".yaml"}:
                    checked_files.append(p)

    for fpath in checked_files:
        text = fpath.read_text(encoding="utf-8")
        assert not _has_cjk(text), f"CJK character detected in {fpath.relative_to(REPO_ROOT)}"
