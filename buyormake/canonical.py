"""Canonical JSON, hashing, atomic download, and RDKit structure normalization helpers."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import httpx
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

_UNCHARGER = rdMolStandardize.Uncharger()


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


def file_md5(path: Path | str, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(Path(path), "rb") as f:
        while block := f.read(chunk):
            h.update(block)
    return h.hexdigest()


def _verify_existing_file(
    path: Path,
    *,
    expected_sha256: Optional[str] = None,
    expected_md5: Optional[str] = None,
    expected_size: Optional[int] = None,
) -> bool:
    if not path.exists():
        return False
    size = path.stat().st_size
    if size <= 0:
        return False
    if expected_size is not None and size != expected_size:
        return False
    if expected_sha256 is not None and file_sha256(path).lower() != expected_sha256.lower():
        return False
    if expected_md5 is not None and file_md5(path).lower() != expected_md5.lower():
        return False
    # Never skip an unverified existing file based solely on arbitrary size unless at least
    # one integrity expectation was checked or caller opted into reuse_unverified.
    return expected_sha256 is not None or expected_md5 is not None or expected_size is not None


def atomic_download(
    url: str,
    dest: Path | str,
    *,
    expected_sha256: Optional[str] = None,
    expected_md5: Optional[str] = None,
    expected_size: Optional[int] = None,
    timeout: float = 120.0,
    force: bool = False,
    client: Optional[httpx.Client] = None,
) -> str:
    """Stream `url` to `<dest>.part`, verify checksums/lengths, and `os.replace` atomically.

    Returns the lowercase hex SHA-256 digest of the final file at `dest`.
    """
    target = Path(dest)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not force and _verify_existing_file(
        target,
        expected_sha256=expected_sha256,
        expected_md5=expected_md5,
        expected_size=expected_size,
    ):
        return file_sha256(target)

    part_path = target.with_suffix(target.suffix + ".part")
    owns_client = client is None
    http = client or httpx.Client(follow_redirects=True, timeout=timeout)
    try:
        with http.stream("GET", url) as resp:
            resp.raise_for_status()
            content_len_hdr = resp.headers.get("Content-Length")
            expected_wire_len = (
                int(content_len_hdr)
                if (content_len_hdr and content_len_hdr.isdigit() and not resp.headers.get("Content-Encoding"))
                else None
            )
            written = 0
            with part_path.open("wb") as fout:
                for chunk in resp.iter_bytes():
                    fout.write(chunk)
                    written += len(chunk)

        if written == 0:
            raise ValueError(f"Downloaded empty payload (0 bytes) from {url}")
        if expected_wire_len is not None and written != expected_wire_len:
            raise ValueError(
                f"Incomplete download from {url}: got {written} bytes, Content-Length={expected_wire_len}"
            )
        if expected_size is not None and written != expected_size:
            raise ValueError(
                f"Size mismatch for {target.name}: got {written} bytes, expected {expected_size}"
            )

        actual_sha256 = file_sha256(part_path)
        if expected_sha256 is not None and actual_sha256.lower() != expected_sha256.lower():
            raise ValueError(
                f"SHA-256 mismatch for {target.name}: got {actual_sha256}, expected {expected_sha256}"
            )
        if expected_md5 is not None:
            actual_md5 = file_md5(part_path)
            if actual_md5.lower() != expected_md5.lower():
                raise ValueError(
                    f"MD5 mismatch for {target.name}: got {actual_md5}, expected {expected_md5}"
                )

        os.replace(part_path, target)
        return actual_sha256
    except Exception:
        if part_path.exists():
            try:
                part_path.unlink()
            except OSError:
                pass
        raise
    finally:
        if owns_client:
            http.close()


@dataclass(frozen=True)
class NormalizedMol:
    canonical_smiles: str
    inchikey: str
    formula: str
    mol_weight: float
    num_heavy_atoms: int = 0
    has_permanent_charge: bool = False
    mol: Optional[Chem.Mol] = field(default=None, repr=False, compare=False, hash=False)

    @property
    def connectivity_block(self) -> str:
        """14-character IUPAC InChIKey skeletal connectivity block."""
        return self.inchikey[:14] if len(self.inchikey) >= 14 else self.inchikey


def normalize_smiles(smiles: str, *, desalt: bool = False) -> Optional[NormalizedMol]:
    """Parse and canonicalize SMILES via RDKit. Returns None for invalid/empty input.

    When `desalt=True`, strips pharmaceutical acid/base salts, hydrates, and counterions via
    `rdMolStandardize.FragmentParent` followed by `Uncharger` so salt forms match neutral
    parent catalog entries. Permanent ions whose formal charge cannot be neutralized by proton
    transfer (such as quaternary ammonium or pyridinium salts) retain their counterions so
    electroneutrality and counterion stoichiometry/molecular weight are preserved.
    """
    if not isinstance(smiles, str):
        return None
    text = smiles.strip()
    if not text:
        return None
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        return None
    if any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms()):
        mol = Chem.Mol(mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)

    has_perm_charge = False
    has_formal_charge = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())

    if desalt:
        try:
            orig_charge = Chem.GetFormalCharge(mol)
            parent = rdMolStandardize.FragmentParent(mol)
            if parent is not None and parent.GetNumAtoms() > 0:
                uncharged = _UNCHARGER.uncharge(parent)
                cand = uncharged if uncharged is not None else parent
                cand_charge = Chem.GetFormalCharge(cand)
                if cand_charge == 0:
                    mol = cand
                else:
                    has_perm_charge = True
                    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
                    if len(frags) <= 1 or abs(cand_charge) <= abs(orig_charge):
                        mol = cand
                    else:
                        # Strip neutral solvates/hydrates while retaining counterions that balance permanent charge
                        charged_frags = [f for f in frags if Chem.GetFormalCharge(f) != 0]
                        if (
                            len(charged_frags) >= 2
                            and len(charged_frags) < len(frags)
                            and sum(Chem.GetFormalCharge(f) for f in charged_frags) == orig_charge
                        ):
                            combined = charged_frags[0]
                            for extra in charged_frags[1:]:
                                combined = Chem.CombineMols(combined, extra)
                            mol = combined
        except Exception:
            pass
    elif has_formal_charge:
        try:
            parent = rdMolStandardize.FragmentParent(mol)
            if parent is not None and parent.GetNumAtoms() > 0:
                uncharged = _UNCHARGER.uncharge(parent)
                cand = uncharged if uncharged is not None else parent
                has_perm_charge = Chem.GetFormalCharge(cand) != 0
        except Exception:
            pass

    can = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    inchikey = Chem.MolToInchiKey(mol)
    if not inchikey:
        return None
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = float(Descriptors.MolWt(mol))
    n_heavy = int(mol.GetNumHeavyAtoms())
    return NormalizedMol(
        canonical_smiles=can,
        inchikey=inchikey,
        formula=formula,
        mol_weight=mw,
        num_heavy_atoms=n_heavy,
        has_permanent_charge=has_perm_charge,
        mol=mol,
    )
