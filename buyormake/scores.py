"""Synthetic accessibility scorers: SAscore (RDKit Contrib), SCScore (MIT NumPy), optional RAscore."""

from __future__ import annotations

import gzip
import json
import math
import os
import sys
from pathlib import Path
from typing import Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, RDConfig

from .schema import SAScores

DEFAULT_SCSCORE_WEIGHTS = Path("data/scscore/model.ckpt-10654.as_numpy.json.gz")

_sascorer_mod = None
_morgan_gen_cache: dict[tuple[int, int], object] = {}


def _get_sascorer():
    global _sascorer_mod
    if _sascorer_mod is None:
        sa_dir = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_dir not in sys.path:
            sys.path.insert(0, sa_dir)
        import sascorer  # type: ignore

        _sascorer_mod = sascorer
    return _sascorer_mod


def compute_sascore(smiles: str, *, mol: Optional[Chem.Mol] = None) -> Optional[float]:
    """Compute Ertl-Schuffenhauer SAscore (1=easy, 10=hard) via RDKit Contrib (BSD-3-Clause)."""
    parsed = mol if mol is not None else Chem.MolFromSmiles(smiles)
    if parsed is None:
        return None
    scorer = _get_sascorer()
    return round(float(scorer.calculateScore(parsed)), 4)


class SCScoreModel:
    """Standalone NumPy implementation of Coley et al. SCScore (MIT License)."""

    def __init__(self, weight_path: Path | str = DEFAULT_SCSCORE_WEIGHTS) -> None:
        self.weight_path = Path(weight_path)
        self.fp_len = 1024
        self.fp_rad = 2
        self.score_scale = 5.0
        self.vars: list[np.ndarray] = []
        self._loaded = False

    def load(self) -> bool:
        if self._loaded:
            return True
        if not self.weight_path.exists():
            return False
        with gzip.open(self.weight_path, "rt", encoding="utf-8") as fin:
            raw = json.load(fin)
        self.vars = [np.array(x, dtype=np.float32) for x in raw]
        self._loaded = True
        return True

    def _mol_to_fp(self, mol: Chem.Mol) -> np.ndarray:
        try:
            from rdkit.Chem import rdFingerprintGenerator

            key = (self.fp_rad, self.fp_len)
            gen = _morgan_gen_cache.get(key)
            if gen is None:
                gen = rdFingerprintGenerator.GetMorganGenerator(
                    radius=self.fp_rad, fpSize=self.fp_len, includeChirality=True
                )
                _morgan_gen_cache[key] = gen
            return gen.GetFingerprintAsNumPy(mol).astype(np.float32)  # type: ignore[attr-defined]
        except Exception:
            bv = AllChem.GetMorganFingerprintAsBitVect(
                mol, self.fp_rad, nBits=self.fp_len, useChirality=True
            )
            return np.array(bv, dtype=np.float32)

    def score(self, smiles: str, *, mol: Optional[Chem.Mol] = None) -> Optional[float]:
        if not self.load():
            return None
        parsed = mol if mol is not None else Chem.MolFromSmiles(smiles)
        if parsed is None:
            return None
        x = self._mol_to_fp(parsed)
        if float(x.sum()) == 0.0:
            # Minimum complexity floor in SCScore's [1.0, 5.0] domain (e.g. [H][H])
            return 1.0
        for i in range(0, len(self.vars), 2):
            W = self.vars[i]
            b = self.vars[i + 1]
            x = np.matmul(x, W) + b
            if i < len(self.vars) - 2:
                x = np.maximum(x, 0.0)
        scalar = float(np.asarray(x).item())
        if not math.isfinite(scalar):
            return None
        # Overflow-safe logistic sigmoid: 1 / (1 + exp(-z)) == 0.5 * (1 + tanh(0.5 * z))
        sigmoid = 0.5 * (1.0 + math.tanh(0.5 * scalar))
        val = 1.0 + (self.score_scale - 1.0) * sigmoid
        return round(val, 4)


_scscore_models: dict[Path, SCScoreModel] = {}


def compute_scscore(
    smiles: str,
    weight_path: Path | str = DEFAULT_SCSCORE_WEIGHTS,
    *,
    mol: Optional[Chem.Mol] = None,
) -> Optional[float]:
    """Compute Coley SCScore (1..5) if weights exist locally; otherwise return None."""
    p = Path(weight_path)
    key = p.resolve() if p.exists() else p
    model = _scscore_models.get(key)
    if model is None:
        model = SCScoreModel(p)
        if p.exists():
            _scscore_models[key] = model
    return model.score(smiles, mol=mol)


def compute_scores(
    smiles: str,
    scscore_weights: Path | str = DEFAULT_SCSCORE_WEIGHTS,
    *,
    mol: Optional[Chem.Mol] = None,
) -> SAScores:
    """Compute available SA scores; missing optional models remain None."""
    parsed = mol if mol is not None else Chem.MolFromSmiles(smiles)
    if parsed is None:
        return SAScores(sascore=None, scscore=None, rascore=None)
    return SAScores(
        sascore=compute_sascore(smiles, mol=parsed),
        scscore=compute_scscore(smiles, weight_path=scscore_weights, mol=parsed),
        rascore=None,
    )
