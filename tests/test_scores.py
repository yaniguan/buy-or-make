"""Tests for synthetic accessibility scorers (SAscore, SCScore)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np

from buyormake.scores import SCScoreModel, compute_sascore, compute_scores


def test_sascore_valid_and_invalid() -> None:
    s = compute_sascore("CC(=O)Oc1ccccc1C(=O)O")
    assert s is not None
    assert 1.0 <= s <= 10.0
    assert compute_sascore("not_valid") is None


def test_scscore_missing_weights_returns_none(tmp_path: Path) -> None:
    scores = compute_scores("CCO", scscore_weights=tmp_path / "missing.json.gz")
    assert scores.sascore is not None
    assert scores.scscore is None
    assert scores.rascore is None


def test_scscore_loads_json_gz(tmp_path: Path) -> None:
    # Minimal 1-layer weights matching 1024-bit FP -> scalar
    w = np.zeros((1024, 1), dtype=np.float32).tolist()
    b = [0.0]
    p = tmp_path / "tiny_scscore.json.gz"
    with gzip.open(p, "wt", encoding="utf-8") as f:
        json.dump([w, b], f)

    model = SCScoreModel(p)
    val = model.score("CC(=O)Oc1ccccc1C(=O)O")
    assert val is not None
    assert 1.0 <= val <= 5.0
