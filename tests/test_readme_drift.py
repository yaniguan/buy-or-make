"""Regression guard verifying zero drift between `results/retrosynthesis_benchmark.json` and `README.md`,
plus repository-wide pure-English (zero CJK) and gitignore hygiene checks.
"""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

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


def test_retrosynthesis_benchmark_json_schema_and_metrics() -> None:
    assert BENCHMARK_JSON.exists(), (
        "Missing committed results/retrosynthesis_benchmark.json; run `python scripts/benchmark_retrosynthesis.py --update-readme`"
    )
    report = json.loads(BENCHMARK_JSON.read_text(encoding="utf-8"))
    assert report["benchmark_schema"] == "buy-or-make/retrosynthesis-benchmark/v1"
    assert report["data_source"] == "open_stopping_stock"

    # Verify exact 64-hex SHA-256 dataset digests for PaRoutes n1/n5, USPTO-50K, SynPlanner, and SCScore
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
    overall = ms["overall"]
    for m_key in (
        "stock_only",
        "lookup_1step_smarts_2step",
        "retro_star_max4",
        "retro_star_max5",
    ):
        assert m_key in overall
        assert overall[m_key]["n_targets"] == ms["n_targets_total"]
        assert 0.0 <= overall[m_key]["solve_rate"] <= 1.0
        assert overall[m_key]["median_time_ms"] >= 0.0
        assert overall[m_key]["solved_leaf_stock_completeness"] == 1.0

    # Retro* must outperform 2-step greedy SMARTS and stock-only lookup on real multi-step targets
    assert overall["retro_star_max5"]["solve_rate"] >= overall["retro_star_max4"]["solve_rate"]
    assert overall["retro_star_max5"]["solve_rate"] > overall["lookup_1step_smarts_2step"]["solve_rate"]
    assert overall["lookup_1step_smarts_2step"]["solve_rate"] > overall["stock_only"]["solve_rate"]

    # Held-out USPTO-50K single-step top-k accuracy monotonicity
    u50 = report["uspto_50k_single_step"]
    assert u50["n_test_targets"] >= 200
    ext_m = u50["methods"]["retro_star_lookup_trainval_plus_44_smarts"]
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

    scan_dirs = ("buyormake", "community_extensions", "docs", "examples", "lean", "results", "schema", "scripts", "tests")
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
