"""Canonical JSON, hashing, and RDKit structure normalization helpers."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.error")


def canonical_json(obj) -> bytes:
    """Sorted keys, no whitespace, UTF-8. Floats use Python's shortest repr."""
    return json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_sha256(path: Path | str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(Path(path), "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


@dataclass(frozen=True)
class NormalizedMol:
    canonical_smiles: str
    inchikey: str
    formula: str
    mol_weight: float


def normalize_smiles(smiles: str) -> Optional[NormalizedMol]:
    """Parse and canonicalize SMILES via RDKit. Returns None for invalid/empty input."""
    if not isinstance(smiles, str):
        return None
    text = smiles.strip()
    if not text:
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    can = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    inchikey = Chem.MolToInchiKey(mol)
    if not inchikey:
        return None
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = float(Descriptors.MolWt(mol))
    return NormalizedMol(
        canonical_smiles=can,
        inchikey=inchikey,
        formula=formula,
        mol_weight=mw,
    )
