"""Tests for synthetic accessibility scorers (SAscore, SCScore)."""

from __future__ import annotations

import gzip
import json
import math
from pathlib import Path

import numpy as np
import pytest

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
    # 2-layer net whose output does not depend on the fingerprint, so the expected value can be
    # computed by hand:
    #   hidden = relu(x @ 0 + [1.0, -1.0]) = [1.0, 0.0]
    #   scalar = hidden @ [[-2.0], [3.0]] + 0.5 = -1.5
    #   score  = 1 + 4 / (1 + exp(1.5)) = 1.7297
    # Dropping the hidden ReLU gives -4.5 -> 1.0439; applying ReLU to the output gives 0.0 -> 3.0.
    w1 = np.zeros((1024, 2), dtype=np.float32).tolist()
    b1 = [1.0, -1.0]
    w2 = [[-2.0], [3.0]]
    b2 = [0.5]
    p = tmp_path / "tiny_scscore.json.gz"
    with gzip.open(p, "wt", encoding="utf-8") as f:
        json.dump([w1, b1, w2, b2], f)

    model = SCScoreModel(p)
    val = model.score("CC(=O)Oc1ccccc1C(=O)O")
    assert val is not None
    expected = 1.0 + 4.0 / (1.0 + math.exp(1.5))
    assert val == pytest.approx(expected, abs=1e-4)
    assert val == pytest.approx(1.7297, abs=1e-4)
