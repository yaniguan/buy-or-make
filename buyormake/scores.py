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


def _get_sascorer():
    global _sascorer_mod
    if _sascorer_mod is None:
        sa_dir = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_dir not in sys.path:
            sys.path.insert(0, sa_dir)
        import sascorer  # type: ignore

        _sascorer_mod = sascorer
    return _sascorer_mod


def compute_sascore(smiles: str) -> Optional[float]:
    """Compute Ertl-Schuffenhauer SAscore (1=easy, 10=hard) via RDKit Contrib (BSD-3-Clause)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scorer = _get_sascorer()
    return round(float(scorer.calculateScore(mol)), 4)


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

            gen = rdFingerprintGenerator.GetMorganGenerator(
                radius=self.fp_rad, fpSize=self.fp_len, includeChirality=True
            )
            return gen.GetFingerprintAsNumPy(mol).astype(np.float32)
        except Exception:
            bv = AllChem.GetMorganFingerprintAsBitVect(
                mol, self.fp_rad, nBits=self.fp_len, useChirality=True
            )
            return np.array(bv, dtype=np.float32)

    def score(self, smiles: str) -> Optional[float]:
        if not self.load():
            return None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        x = self._mol_to_fp(mol)
        if float(x.sum()) == 0.0:
            return 0.0
        for i in range(0, len(self.vars), 2):
            W = self.vars[i]
            b = self.vars[i + 1]
            x = np.matmul(x, W) + b
            if i < len(self.vars) - 2:
                x = np.maximum(x, 0.0)
        scalar = float(np.asarray(x).item())
        val = 1.0 + (self.score_scale - 1.0) / (1.0 + math.exp(-scalar))
        return round(val, 4)


_default_scscore_model: Optional[SCScoreModel] = None


def compute_scscore(
    smiles: str,
    weight_path: Path | str = DEFAULT_SCSCORE_WEIGHTS,
) -> Optional[float]:
    """Compute Coley SCScore (1..5) if weights exist locally; otherwise return None."""
    global _default_scscore_model
    p = Path(weight_path)
    if Path(weight_path) == DEFAULT_SCSCORE_WEIGHTS:
        if _default_scscore_model is None:
            _default_scscore_model = SCScoreModel(p)
        return _default_scscore_model.score(smiles)
    return SCScoreModel(p).score(smiles)


def compute_scores(
    smiles: str,
    scscore_weights: Path | str = DEFAULT_SCSCORE_WEIGHTS,
) -> SAScores:
    """Compute available SA scores; missing optional models remain None."""
    return SAScores(
        sascore=compute_sascore(smiles),
        scscore=compute_scscore(smiles, weight_path=scscore_weights),
        rascore=None,
    )
