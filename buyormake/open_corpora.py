"""On-demand fetchers and multi-format importers for 2024-2026 open chemistry corpora.

All downloaded dumps land in gitignored `data/open_corpora/<source_id>/` via checksum-verified
atomic `.part` downloads (`buyormake.canonical.atomic_download`). Stopping-stock datasets index
into `data/stock.sqlite`, reaction/route benchmarks index into `data/reactions.sqlite`, and
non-vendor reference compound corpora (`ChEMBL`, `SureChEMBL`, `COCONUT 2.0`) import into
`data/reference_compounds.sqlite` (or sealed `AssayfoldRecord` envelopes) so unpriced reference
molecules are never silently conflated with commercial vendor stock.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import sqlite3
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Literal, Optional, TextIO

import httpx

from .canonical import NormalizedMol, atomic_download, normalize_smiles
from .cost import EST_COST_METHOD
from .routes import parse_reaction_smiles_for_indexing
from .schema import (
    AssayfoldRecord,
    Molecule,
    Provenance,
    Route,
    RouteStep,
    make_assayfold_record_id,
    utc_now,
)
from .sources import SOURCE_BY_ID

DEFAULT_OPEN_CORPORA_DIR = Path("data/open_corpora")
DEFAULT_STOCK_DB = Path("data/stock.sqlite")
DEFAULT_REFERENCE_COMPOUNDS_DB = Path("data/reference_compounds.sqlite")
DEFAULT_REACTION_DB = Path("data/reactions.sqlite")

CorpusKind = Literal["stock", "reference_compounds", "reactions", "routes", "targets"]


@dataclass(frozen=True)
class CorpusArtifact:
    artifact_key: str
    source_id: str
    description: str
    url: str
    filename: str
    kind: CorpusKind
    expected_sha256: Optional[str] = None
    expected_md5: Optional[str] = None
    expected_size: Optional[int] = None


CORPUS_ARTIFACTS: dict[str, CorpusArtifact] = {
    "paroutes_stock_n1": CorpusArtifact(
        artifact_key="paroutes_stock_n1",
        source_id="paroutes",
        description="PaRoutes 2.0 n1 stopping-stock building block SMILES (CC-BY-4.0)",
        url="https://zenodo.org/api/records/7341155/files/stock_n1.txt/content",
        filename="stock_n1.txt",
        kind="stock",
        expected_md5="cfa7e97b0398c7e88e46c02ae5876edf",
        expected_size=368_802,
    ),
    "paroutes_stock_n5": CorpusArtifact(
        artifact_key="paroutes_stock_n5",
        source_id="paroutes",
        description="PaRoutes 2.0 n5 stopping-stock building block SMILES (CC-BY-4.0)",
        url="https://zenodo.org/api/records/7341155/files/stock_n5.txt/content",
        filename="stock_n5.txt",
        kind="stock",
        expected_md5="741a804816fd0484a683b58846e142bb",
        expected_size=360_970,
    ),
    "paroutes_routes_n1": CorpusArtifact(
        artifact_key="paroutes_routes_n1",
        source_id="paroutes",
        description="PaRoutes 2.0 n1 10,000 multi-step reference route trees (CC-BY-4.0)",
        url="https://zenodo.org/api/records/7341155/files/ref_routes_n1.json/content",
        filename="ref_routes_n1.json",
        kind="routes",
        expected_md5="8f11ca81d6d184539fa4d2c5bde215af",
        expected_size=64_187_987,
    ),
    "paroutes_routes_n5": CorpusArtifact(
        artifact_key="paroutes_routes_n5",
        source_id="paroutes",
        description="PaRoutes 2.0 n5 10,000 multi-step reference route trees (CC-BY-4.0)",
        url="https://zenodo.org/api/records/7341155/files/ref_routes_n5.json/content",
        filename="ref_routes_n5.json",
        kind="routes",
        expected_md5="1da36cbada13175668c66546db6fd952",
        expected_size=82_120_928,
    ),
    "paroutes_reactions": CorpusArtifact(
        artifact_key="paroutes_reactions",
        source_id="paroutes",
        description="PaRoutes 2.0 curated USPTO reaction SMILES CSV (CC-BY-4.0)",
        url="https://zenodo.org/api/records/7341155/files/selected_reactions.csv/content",
        filename="selected_reactions.csv",
        kind="reactions",
        expected_md5="0e77ccbfd0ff650afe09148d13cb04ee",
        expected_size=582_648_364,
    ),
    "chembl_37_chemreps": CorpusArtifact(
        artifact_key="chembl_37_chemreps",
        source_id="chembl",
        description="ChEMBL 37 (May 2026) canonical SMILES & InChIKey TSV (CC-BY-SA-3.0)",
        url="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_chemreps.txt.gz",
        filename="chembl_37_chemreps.txt.gz",
        kind="reference_compounds",
        expected_sha256="ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7",
        expected_size=292_540_141,
    ),
    "surechembl_compounds_parquet": CorpusArtifact(
        artifact_key="surechembl_compounds_parquet",
        source_id="surechembl",
        description="SureChEMBL bulk_data/latest patent compounds Parquet (CC-BY-4.0)",
        url="https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/latest/compounds.parquet",
        filename="compounds.parquet",
        kind="reference_compounds",
        expected_size=4_281_834_894,
    ),
    "coconut_csv_lite_2026_11": CorpusArtifact(
        artifact_key="coconut_csv_lite_2026_11",
        source_id="coconut",
        description="COCONUT 2.0 monthly natural products lite CSV archive (2026-11, CC0-1.0)",
        url="https://coconut.s3.uni-jena.de/prod/downloads/2026-11/coconut_csv_lite-11-2026.zip",
        filename="coconut_csv_lite-11-2026.zip",
        kind="reference_compounds",
        expected_size=199_108_039,
    ),
    "coconut_2_0_zenodo_csv": CorpusArtifact(
        artifact_key="coconut_2_0_zenodo_csv",
        source_id="coconut",
        description="COCONUT 2.0 archived Zenodo 13897048 CSV snapshot",
        url="https://zenodo.org/api/records/13897048/files/coconut-10-2024.csv.zip/content",
        filename="coconut-10-2024.csv.zip",
        kind="reference_compounds",
        expected_md5="1d49f92a26f20eeabe449980b8c17307",
        expected_size=188_339_585,
    ),
    "orderly_retro_test": CorpusArtifact(
        artifact_key="orderly_retro_test",
        source_id="orderly",
        description="ORDerly cleaned retrosynthesis test split Parquet (Figshare 23298467, CC-BY-4.0)",
        url="https://ndownloader.figshare.com/files/44413049",
        filename="orderly_retro_test.parquet",
        kind="reactions",
        expected_md5="3dc796d5fbd61de00baf9945e13d4fb1",
        expected_size=53_303_576,
    ),
    "orderly_retro_train": CorpusArtifact(
        artifact_key="orderly_retro_train",
        source_id="orderly",
        description="ORDerly cleaned retrosynthesis train split Parquet (Figshare 23298467, CC-BY-4.0)",
        url="https://ndownloader.figshare.com/files/44413061",
        filename="orderly_retro_train.parquet",
        kind="reactions",
        expected_md5="a97ae485f100af2903069dab85cc4b07",
        expected_size=520_416_959,
    ),
    "syntheseus_uspto_50k": CorpusArtifact(
        artifact_key="syntheseus_uspto_50k",
        source_id="syntheseus_uspto_50k",
        description="Syntheseus USPTO-50K single-step retrosynthesis benchmark (Figshare 25459573, MIT)",
        url="https://ndownloader.figshare.com/files/45206101",
        filename="uspto50k.zip",
        kind="reactions",
        expected_md5="b4c1128ed1bf663f48ff17a1bc2e6517",
        expected_size=19_940_280,
    ),
    "syntheseus_retro_star_inventory": CorpusArtifact(
        artifact_key="syntheseus_retro_star_inventory",
        source_id="syntheseus_uspto_50k",
        description="Syntheseus Retro* building-block stock inventory CSV.GZ (Figshare 25376728, MIT)",
        url="https://ndownloader.figshare.com/files/44956075",
        filename="ORIGINAL-inventory.csv.gz",
        kind="stock",
        expected_md5="48b61a98b9d62876c9d082d2014c6e7a",
        expected_size=273_194_394,
    ),
    "synplanner_building_blocks": CorpusArtifact(
        artifact_key="synplanner_building_blocks",
        source_id="synplanner_data",
        description="SynPlanner building-block SMILES on Hugging Face Hub (MIT)",
        url="https://huggingface.co/Laboratoire-De-Chemoinformatique/SynPlanner/resolve/main/building_blocks/building_blocks_em_sa_ln.smi",
        filename="building_blocks_em_sa_ln.smi",
        kind="stock",
        expected_size=5_784_208,
    ),
    "synplanner_targets_sa_1_5_2_5": CorpusArtifact(
        artifact_key="synplanner_targets_sa_1_5_2_5",
        source_id="synplanner_data",
        description="SynPlanner SAscore [1.5, 2.5) retrosynthesis benchmark targets (MIT)",
        url="https://huggingface.co/Laboratoire-De-Chemoinformatique/SynPlanner/resolve/main/benchmarks/sascore/targets_with_sascore_1.5_2.5.smi",
        filename="targets_with_sascore_1.5_2.5.smi",
        kind="targets",
        expected_size=4_361,
    ),
    "synplanner_uspto_standardized": CorpusArtifact(
        artifact_key="synplanner_uspto_standardized",
        source_id="synplanner_data",
        description="SynPlanner standardized USPTO reaction corpus ZIP on Hugging Face Hub (MIT)",
        url="https://huggingface.co/Laboratoire-De-Chemoinformatique/SynPlanner/resolve/main/uspto/uspto_standardized.zip",
        filename="uspto_standardized.zip",
        kind="reactions",
        expected_sha256="c40af53513925e3d47a547ad6fe6bc475e2ae07ec16f43f540e6cdd8bd3b9908",
        expected_size=118_758_231,
    ),
}


def download_corpus_artifact(
    artifact_key: str,
    out_root: Path | str = DEFAULT_OPEN_CORPORA_DIR,
    *,
    expected_sha256: Optional[str] = None,
    force: bool = False,
    timeout: float = 180.0,
    client: Optional[httpx.Client] = None,
) -> Path:
    """Download a registered open-corpus artifact into gitignored `data/open_corpora/<source_id>/`."""
    if artifact_key not in CORPUS_ARTIFACTS:
        raise KeyError(
            f"Unknown corpus artifact {artifact_key!r}; expected one of {sorted(CORPUS_ARTIFACTS)}"
        )
    spec = CORPUS_ARTIFACTS[artifact_key]
    if spec.source_id not in SOURCE_BY_ID:
        raise KeyError(f"Artifact {artifact_key!r} references unregistered source {spec.source_id!r}")
    safe_name = Path(spec.filename).name
    if safe_name != spec.filename or ".." in spec.filename:
        raise ValueError(f"Unsafe artifact filename {spec.filename!r}")
    dest = Path(out_root) / spec.source_id / safe_name
    atomic_download(
        spec.url,
        dest,
        expected_sha256=expected_sha256 or spec.expected_sha256,
        expected_md5=spec.expected_md5,
        expected_size=spec.expected_size,
        timeout=timeout,
        force=force,
        client=client,
    )
    return dest


# ---------- SQLite schema helpers ----------


def _init_stock_db(conn: sqlite3.Connection) -> None:
    conn.execute("PRAGMA busy_timeout = 30000")
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS stock (
            inchikey TEXT NOT NULL,
            parent_inchikey TEXT,
            smiles TEXT,
            vendor TEXT NOT NULL,
            catalog_id TEXT,
            pack_value REAL,
            pack_unit TEXT,
            purity TEXT,
            price REAL,
            currency TEXT,
            price_basis TEXT NOT NULL DEFAULT 'unknown',
            lead_time_days INTEGER,
            in_stock INTEGER,
            source_id TEXT NOT NULL,
            source_url TEXT NOT NULL,
            fetched_utc TEXT NOT NULL,
            license TEXT NOT NULL
        )
        """
    )
    cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(stock)").fetchall()}
    if "parent_inchikey" not in cols:
        conn.execute("ALTER TABLE stock ADD COLUMN parent_inchikey TEXT")
    if "smiles" not in cols:
        conn.execute("ALTER TABLE stock ADD COLUMN smiles TEXT")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_stock_ik ON stock(inchikey)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_stock_parent_ik ON stock(parent_inchikey)")
    conn.commit()


def _init_reaction_db(conn: sqlite3.Connection) -> None:
    conn.execute("PRAGMA busy_timeout = 30000")
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS reactions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            product_inchikey TEXT NOT NULL,
            product_smiles TEXT NOT NULL,
            reactants_json TEXT NOT NULL,
            agents_json TEXT NOT NULL DEFAULT '[]',
            reaction_smiles TEXT NOT NULL,
            reference TEXT NOT NULL,
            source TEXT NOT NULL
        )
        """
    )
    cols = {str(r[1]) for r in conn.execute("PRAGMA table_info(reactions)").fetchall()}
    if "agents_json" not in cols:
        conn.execute("ALTER TABLE reactions ADD COLUMN agents_json TEXT NOT NULL DEFAULT '[]'")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_rxn_prod_ik ON reactions(product_inchikey)")
    conn.commit()


def _open_text_streams(path: Path, suffixes: tuple[str, ...] = ()) -> Iterator[tuple[str, TextIO]]:
    """Yield `(member_name, text_stream)` from plain text, `.gz`, or `.zip` files."""
    lower = path.name.lower()
    if lower.endswith(".zip"):
        with zipfile.ZipFile(path, "r") as zf:
            for info in zf.infolist():
                if info.is_dir():
                    continue
                inner_name = Path(info.filename).name
                if inner_name.startswith("."):
                    continue
                if suffixes and not inner_name.lower().endswith(suffixes):
                    continue
                with zf.open(info, "r") as raw:
                    if inner_name.lower().endswith(".gz"):
                        with gzip.open(raw, "rt", encoding="utf-8", errors="replace") as gz_stream:
                            yield inner_name, gz_stream
                    else:
                        wrapper = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                        yield inner_name, wrapper
        return

    if lower.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as gz_stream:
            yield path.name, gz_stream
        return

    with path.open("r", encoding="utf-8", errors="replace") as f:
        yield path.name, f


def _parent_ik_if_salted(nm: NormalizedMol) -> Optional[str]:
    if "." in nm.canonical_smiles or "+" in nm.canonical_smiles or "-" in nm.canonical_smiles:
        p_nm = normalize_smiles(nm.canonical_smiles, desalt=True)
        if p_nm is not None and p_nm.inchikey != nm.inchikey:
            return p_nm.inchikey
    return None


def _insert_unpriced_molecules(
    entries: Iterable[tuple[NormalizedMol, Optional[str]]],
    *,
    db_path: Path,
    vendor: str,
    source_id: str,
    source_url: str,
    license_spdx: str,
    in_stock: Optional[bool],
    max_rows: int = 0,
) -> int:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path), timeout=30.0)
    _init_stock_db(conn)
    ts = utc_now()
    in_stock_int = None if in_stock is None else (1 if in_stock else 0)
    sql = (
        "INSERT INTO stock (inchikey, parent_inchikey, smiles, vendor, catalog_id, "
        "price_basis, in_stock, source_id, source_url, fetched_utc, license) "
        "VALUES (?, ?, ?, ?, ?, 'unknown', ?, ?, ?, ?, ?)"
    )
    batch: list[tuple[Any, ...]] = []
    inserted = 0
    try:
        for nm, cat_id in entries:
            parent_ik = _parent_ik_if_salted(nm)
            batch.append(
                (
                    nm.inchikey,
                    parent_ik,
                    nm.canonical_smiles,
                    vendor,
                    cat_id,
                    in_stock_int,
                    source_id,
                    source_url,
                    ts,
                    license_spdx,
                )
            )
            inserted += 1
            if len(batch) >= 2000:
                conn.executemany(sql, batch)
                conn.commit()
                batch.clear()
            if max_rows and inserted >= max_rows:
                break
        if batch:
            conn.executemany(sql, batch)
            conn.commit()
    finally:
        conn.close()
    return inserted


def _index_reaction_pairs(
    pairs: Iterable[tuple[Optional[str], Optional[str]]],
    *,
    db_path: Path,
    step_source: str,
    default_prefix: str,
    max_rows: int = 0,
) -> int:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path), timeout=30.0)
    _init_reaction_db(conn)
    sql = (
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, "
        "agents_json, reaction_smiles, reference, source) VALUES (?, ?, ?, ?, ?, ?, ?)"
    )
    batch: list[tuple[str, str, str, str, str, str, str]] = []
    inserted = 0
    try:
        for idx, (raw_rxn, ref) in enumerate(pairs):
            if not raw_rxn:
                continue
            rxn_smi = str(raw_rxn).strip()
            # Normalize three-part `reactants>agents>products` into `reactants>agents>>products`
            # so `parse_reaction_smiles_for_indexing` routes `agents` to `agents_json`.
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
                if max_rows and inserted >= max_rows:
                    break
            if max_rows and inserted >= max_rows:
                break
        if batch:
            conn.executemany(sql, batch)
            conn.commit()
    finally:
        conn.close()
    return inserted


def _read_parquet_columns(parquet_path: Path) -> dict[str, list[Any]]:
    try:
        import pyarrow.parquet as pq  # type: ignore

        table = pq.read_table(parquet_path)
        return {col: table[col].to_pylist() for col in table.column_names}
    except ImportError:
        try:
            import pandas as pd  # type: ignore

            df = pd.read_parquet(parquet_path)
            return {str(col): df[col].tolist() for col in df.columns}
        except ImportError as exc:
            raise RuntimeError(
                f"Reading Parquet file {parquet_path} requires pyarrow or pandas (`pip install -e .[data]`)."
            ) from exc


def _pick_col(headers: Iterable[str], candidates: tuple[str, ...]) -> Optional[str]:
    by_lower = {h.strip().lower(): h for h in headers if h}
    for cand in candidates:
        if cand.lower() in by_lower:
            return by_lower[cand.lower()]
    return None


# ---------- 1. PaRoutes 2.0 (stock + multi-step route trees + reactions) ----------


def import_paroutes_stock(
    path: Path | str,
    db_path: Path | str = DEFAULT_STOCK_DB,
    *,
    subset: str = "n1",
    max_rows: int = 0,
) -> int:
    """Import PaRoutes `stock_n1.txt` / `stock_n5.txt` SMILES into `stock.sqlite` (`CC-BY-4.0`)."""
    src = SOURCE_BY_ID["paroutes"]
    p = Path(path)

    def _iter() -> Iterator[tuple[NormalizedMol, Optional[str]]]:
        idx = 0
        for _, stream in _open_text_streams(p, suffixes=(".txt", ".smi", ".csv", ".tsv")):
            for raw_line in stream:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if idx == 0 and parts[0].lower() in {"smiles", "inchi_key", "inchikey"}:
                    continue
                nm = normalize_smiles(parts[0], desalt=False)
                if nm is None:
                    continue
                idx += 1
                cat_id = parts[1] if len(parts) > 1 else f"paroutes-{subset}:{idx}"
                yield nm, cat_id

    return _insert_unpriced_molecules(
        _iter(),
        db_path=Path(db_path),
        vendor=f"PaRoutes-{subset}",
        source_id=src.id,
        source_url="https://zenodo.org/records/7341155",
        license_spdx=src.spdx_license or src.license,
        in_stock=True,
        max_rows=max_rows,
    )


def _walk_paroutes_node(
    node: dict[str, Any],
    steps: list[RouteStep],
    leaves: list[str],
    leaves_ok: list[bool],
) -> None:
    node_type = node.get("type", "mol")
    if node_type == "mol":
        raw_smi = str(node.get("smiles") or "").strip()
        nm = normalize_smiles(raw_smi, desalt=False) if raw_smi else None
        canon_smi = nm.canonical_smiles if nm is not None else raw_smi
        children = [c for c in (node.get("children") or []) if isinstance(c, dict)]
        rxn_children = [c for c in children if c.get("type") == "reaction"]
        if not rxn_children and children:
            rxn_children = [children[0]]
        if not rxn_children:
            if canon_smi:
                leaves.append(canon_smi)
                leaves_ok.append(bool(node.get("in_stock", False)))
            return

        best_branch: Optional[tuple[tuple[int, int, int], list[RouteStep], list[str], list[bool]]] = None
        for rxn_node in rxn_children:
            b_steps: list[RouteStep] = []
            b_leaves: list[str] = []
            b_ok: list[bool] = []
            mol_children = [
                c for c in (rxn_node.get("children") or []) if isinstance(c, dict) and c.get("type", "mol") == "mol"
            ]
            if not mol_children:
                continue
            reactant_smiles: list[str] = []
            for mc in mol_children:
                mc_raw = str(mc.get("smiles") or "").strip()
                mc_nm = normalize_smiles(mc_raw, desalt=False) if mc_raw else None
                mc_smi = mc_nm.canonical_smiles if mc_nm is not None else mc_raw
                if mc_smi:
                    reactant_smiles.append(mc_smi)
            if not reactant_smiles or not canon_smi:
                continue
            meta = rxn_node.get("metadata") if isinstance(rxn_node.get("metadata"), dict) else {}
            ref = str(
                meta.get("ID")
                or meta.get("id")
                or meta.get("patent_id")
                or meta.get("template_hash")
                or rxn_node.get("smiles")
                or "paroutes:step"
            )
            rxn_smi = str(
                meta.get("rsmi")
                or meta.get("reaction_smiles")
                or f"{'.'.join(reactant_smiles)}>>{canon_smi}"
            )
            if ">>" not in rxn_smi:
                rxn_smi = f"{'.'.join(reactant_smiles)}>>{canon_smi}"
            b_steps.append(
                RouteStep(
                    product=canon_smi,
                    reactants=reactant_smiles,
                    reaction_smiles=rxn_smi,
                    template_or_reference=ref,
                    source="paroutes_lookup",
                )
            )
            for mc in mol_children:
                _walk_paroutes_node(mc, b_steps, b_leaves, b_ok)
            all_solved = 1 if (b_ok and all(b_ok)) else 0
            rank_key = (-all_solved, len(b_steps), len(b_leaves))
            if best_branch is None or rank_key < best_branch[0]:
                best_branch = (rank_key, b_steps, b_leaves, b_ok)

        if best_branch is not None:
            steps.extend(best_branch[1])
            leaves.extend(best_branch[2])
            leaves_ok.extend(best_branch[3])
        elif canon_smi:
            leaves.append(canon_smi)
            leaves_ok.append(bool(node.get("in_stock", False)))


def parse_paroutes_route_tree(tree: dict[str, Any]) -> Optional[Route]:
    """Convert a PaRoutes JSON route tree dict into a validated `buyormake.schema.Route`."""
    if not isinstance(tree, dict):
        return None
    root = tree.get("dict") if isinstance(tree.get("dict"), dict) else tree
    root_smi = str(root.get("smiles") or "").strip()
    target_nm = normalize_smiles(root_smi, desalt=False)
    if target_nm is None:
        return None
    steps: list[RouteStep] = []
    leaves: list[str] = []
    leaves_ok: list[bool] = []
    _walk_paroutes_node(root, steps, leaves, leaves_ok)
    if not steps or not leaves:
        return None
    return Route(
        target_inchikey=target_nm.inchikey,
        n_steps=len(steps),
        steps=steps,
        leaves=leaves,
        leaves_purchasable=leaves_ok,
        method="paroutes (CC-BY-4.0 reference route)",
        score=None,
        est_cost=None,
        est_cost_method=EST_COST_METHOD,
        est_cost_complete=False,
    )


def load_paroutes_routes(path: Path | str, *, max_routes: int = 0) -> list[Route]:
    """Load multi-step `Route` objects from a PaRoutes `.json` or `.json.gz` file."""
    p = Path(path)
    routes: list[Route] = []
    for _, stream in _open_text_streams(p, suffixes=(".json",)):
        payload = json.load(stream)
        trees = payload if isinstance(payload, list) else [payload]
        for item in trees:
            if not isinstance(item, dict):
                continue
            route = parse_paroutes_route_tree(item)
            if route is not None:
                routes.append(route)
                if max_routes and len(routes) >= max_routes:
                    return routes
    return routes


def index_paroutes_routes(
    path: Path | str,
    db_path: Path | str = DEFAULT_REACTION_DB,
    *,
    max_rows: int = 0,
) -> int:
    """Index reaction steps from PaRoutes reference route trees into `reactions.sqlite`."""
    routes = load_paroutes_routes(path)

    def _pairs() -> Iterator[tuple[Optional[str], Optional[str]]]:
        for r_idx, route in enumerate(routes):
            for s_idx, step in enumerate(route.steps):
                ref = step.template_or_reference or f"paroutes:route{r_idx}:step{s_idx}"
                reactant_set = set(step.reactants)
                extra_agents: list[str] = []
                if step.reaction_smiles and ">>" in step.reaction_smiles:
                    lhs_raw = step.reaction_smiles.split(">>", 1)[0]
                    for part in lhs_raw.replace(">", ".").split("."):
                        part = part.strip()
                        if not part:
                            continue
                        nm = normalize_smiles(part, desalt=False)
                        cand = nm.canonical_smiles if nm is not None else part
                        if cand not in reactant_set:
                            extra_agents.append(cand)
                if extra_agents:
                    rxn = f"{'.'.join(step.reactants)}>{'.'.join(extra_agents)}>>{step.product}"
                else:
                    rxn = f"{'.'.join(step.reactants)}>>{step.product}"
                yield rxn, ref

    return _index_reaction_pairs(
        _pairs(),
        db_path=Path(db_path),
        step_source="paroutes_lookup",
        default_prefix="paroutes",
        max_rows=max_rows,
    )


def index_paroutes_reactions_csv(
    path: Path | str,
    db_path: Path | str = DEFAULT_REACTION_DB,
    *,
    max_rows: int = 0,
) -> int:
    """Index `selected_reactions.csv` (`.csv`, `.tsv`, `.gz`, `.zip`) from PaRoutes into `reactions.sqlite`."""
    p = Path(path)

    def _pairs() -> Iterator[tuple[Optional[str], Optional[str]]]:
        for name, stream in _open_text_streams(p, suffixes=(".csv", ".tsv", ".txt")):
            delim = "\t" if name.lower().endswith(".tsv") else ","
            reader = csv.DictReader(stream, delimiter=delim)
            if not reader.fieldnames:
                continue
            smi_col = _pick_col(
                reader.fieldnames,
                ("rxn_smiles", "reaction_smiles", "rsmi", "reactions", "smiles"),
            )
            id_col = _pick_col(
                reader.fieldnames,
                ("id", "reaction_hash", "pseudo_hash", "patent_number", "patent_id"),
            )
            if smi_col is None:
                continue
            for row in reader:
                yield row.get(smi_col), (row.get(id_col) if id_col else None)

    return _index_reaction_pairs(
        _pairs(),
        db_path=Path(db_path),
        step_source="paroutes_lookup",
        default_prefix="paroutes",
        max_rows=max_rows,
    )


# ---------- 2. ChEMBL 37 (chemreps TSV / GZ -> SQLite & AssayfoldRecord) ----------


def _iter_chembl_rows(path: Path, max_rows: int = 0) -> Iterator[tuple[str, NormalizedMol]]:
    yielded = 0
    for _, stream in _open_text_streams(path, suffixes=(".txt", ".tsv", ".csv")):
        header_cols: Optional[list[str]] = None
        smi_idx = 1
        id_idx = 0
        for raw_line in stream:
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            delim = "\t" if "\t" in line else ","
            parts = line.split(delim)
            if header_cols is None and any(
                "smiles" in c.lower() or "chembl_id" in c.lower() for c in parts
            ):
                header_cols = [c.strip().lower() for c in parts]
                for i, col in enumerate(header_cols):
                    if col in {"canonical_smiles", "smiles"}:
                        smi_idx = i
                    elif col in {"chembl_id", "id", "molecule_chembl_id"}:
                        id_idx = i
                continue
            if len(parts) <= max(smi_idx, id_idx):
                continue
            chembl_id = parts[id_idx].strip()
            raw_smi = parts[smi_idx].strip()
            nm = normalize_smiles(raw_smi, desalt=False)
            if nm is None or not chembl_id:
                continue
            yield chembl_id, nm
            yielded += 1
            if max_rows and yielded >= max_rows:
                return


def import_chembl_chemreps(
    path: Path | str,
    db_path: Path | str = DEFAULT_REFERENCE_COMPOUNDS_DB,
    *,
    release: str = "chembl_37",
    max_rows: int = 0,
) -> int:
    """Import EBI `chembl_XX_chemreps.txt[.gz]` into `reference_compounds.sqlite` (`CC-BY-SA-3.0`)."""
    src = SOURCE_BY_ID["chembl"]
    p = Path(path)
    entries = ((nm, chembl_id) for chembl_id, nm in _iter_chembl_rows(p, max_rows=max_rows))
    return _insert_unpriced_molecules(
        entries,
        db_path=Path(db_path),
        vendor=f"ChEMBL ({release})",
        source_id=src.id,
        source_url="https://www.ebi.ac.uk/chembl/",
        license_spdx=src.spdx_license or src.license,
        in_stock=None,
        max_rows=max_rows,
    )


def parse_chembl_chemreps_records(
    path: Path | str,
    *,
    release: str = "chembl_37",
    max_rows: int = 0,
) -> list[AssayfoldRecord]:
    """Convert ChEMBL `chemreps` rows into sealed `AssayfoldRecord(kind='molecule')` objects."""
    src = SOURCE_BY_ID["chembl"]
    lic = src.spdx_license or src.license
    out: list[AssayfoldRecord] = []
    for chembl_id, nm in _iter_chembl_rows(Path(path), max_rows=max_rows):
        mol = Molecule(
            smiles=nm.canonical_smiles,
            inchikey=nm.inchikey,
            name=chembl_id,
            formula=nm.formula,
        )
        prov = Provenance(
            source_url=f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/",
            used=[src.id],
        )
        rec = AssayfoldRecord(
            record_id=make_assayfold_record_id(src.id, release, chembl_id),
            source_id=src.id,
            dataset_id=release,
            native_id=chembl_id,
            kind="molecule",
            payload=mol,
            license=lic,
            origin="real",
            provenance=prov,
            tags=["corpus:chembl", f"release:{release}", f"license:{lic}"],
        ).sealed()
        out.append(rec)
    return out


# ---------- 3. SureChEMBL (Parquet or TSV/CSV -> reference SQLite & AssayfoldRecord) ----------


def _iter_surechembl_rows(path: Path, max_rows: int = 0) -> Iterator[tuple[str, NormalizedMol]]:
    yielded = 0
    if path.name.lower().endswith(".parquet"):
        cols = _read_parquet_columns(path)
        smi_col = _pick_col(cols.keys(), ("smiles", "canonical_smiles", "structure"))
        id_col = _pick_col(cols.keys(), ("surechembl_id", "schembl_id", "id", "compound_id"))
        if smi_col is None:
            return
        smis = cols[smi_col]
        ids = cols[id_col] if id_col else [f"SCHEMBL_ROW_{i}" for i in range(len(smis))]
        for raw_smi, raw_id in zip(smis, ids):
            if not raw_smi:
                continue
            nm = normalize_smiles(str(raw_smi), desalt=False)
            if nm is None:
                continue
            sc_id = str(raw_id).strip() if raw_id is not None else nm.inchikey
            yield sc_id, nm
            yielded += 1
            if max_rows and yielded >= max_rows:
                return
        return

    for name, stream in _open_text_streams(path, suffixes=(".tsv", ".csv", ".txt")):
        delim = "\t" if name.lower().endswith((".tsv", ".txt")) else ","
        reader = csv.DictReader(stream, delimiter=delim)
        if not reader.fieldnames:
            continue
        smi_col = _pick_col(reader.fieldnames, ("smiles", "canonical_smiles", "structure"))
        id_col = _pick_col(reader.fieldnames, ("surechembl_id", "schembl_id", "id", "compound_id"))
        if smi_col is None:
            continue
        for idx, row in enumerate(reader):
            raw_smi = row.get(smi_col)
            if not raw_smi:
                continue
            nm = normalize_smiles(raw_smi, desalt=False)
            if nm is None:
                continue
            sc_id = (row.get(id_col) or f"SCHEMBL_ROW_{idx}").strip() if id_col else f"SCHEMBL_ROW_{idx}"
            yield sc_id, nm
            yielded += 1
            if max_rows and yielded >= max_rows:
                return


def import_surechembl_compounds(
    path: Path | str,
    db_path: Path | str = DEFAULT_REFERENCE_COMPOUNDS_DB,
    *,
    max_rows: int = 0,
) -> int:
    """Import SureChEMBL Parquet or delimited dump into `reference_compounds.sqlite` (`CC-BY-4.0`)."""
    src = SOURCE_BY_ID["surechembl"]
    entries = ((nm, sc_id) for sc_id, nm in _iter_surechembl_rows(Path(path), max_rows=max_rows))
    return _insert_unpriced_molecules(
        entries,
        db_path=Path(db_path),
        vendor="SureChEMBL",
        source_id=src.id,
        source_url="https://www.surechembl.org/",
        license_spdx=src.spdx_license or src.license,
        in_stock=None,
        max_rows=max_rows,
    )


# ---------- 4. COCONUT 2.0 (CSV / CSV.GZ / ZIP -> reference SQLite & AssayfoldRecord) ----------


def _iter_coconut_rows(
    path: Path, max_rows: int = 0
) -> Iterator[tuple[str, Optional[str], NormalizedMol]]:
    yielded = 0
    for name, stream in _open_text_streams(path, suffixes=(".csv", ".tsv", ".txt")):
        delim = "\t" if name.lower().endswith(".tsv") else ","
        reader = csv.DictReader(stream, delimiter=delim)
        if not reader.fieldnames:
            continue
        smi_col = _pick_col(reader.fieldnames, ("canonical_smiles", "smiles", "sugar_free_smiles"))
        id_col = _pick_col(reader.fieldnames, ("identifier", "coconut_id", "id"))
        name_col = _pick_col(reader.fieldnames, ("name", "iupac_name"))
        if smi_col is None:
            continue
        for idx, row in enumerate(reader):
            raw_smi = row.get(smi_col)
            if not raw_smi:
                continue
            nm = normalize_smiles(raw_smi, desalt=False)
            if nm is None:
                continue
            cid = (row.get(id_col) or f"CNP_ROW_{idx}").strip() if id_col else f"CNP_ROW_{idx}"
            mol_name = (row.get(name_col) or "").strip() or None if name_col else None
            yield cid, mol_name, nm
            yielded += 1
            if max_rows and yielded >= max_rows:
                return


def import_coconut_compounds(
    path: Path | str,
    db_path: Path | str = DEFAULT_REFERENCE_COMPOUNDS_DB,
    *,
    max_rows: int = 0,
) -> int:
    """Import COCONUT 2.0 `.csv`, `.csv.gz`, or `.zip` into `reference_compounds.sqlite` (`CC0-1.0`)."""
    src = SOURCE_BY_ID["coconut"]
    entries = ((nm, cid) for cid, _, nm in _iter_coconut_rows(Path(path), max_rows=max_rows))
    return _insert_unpriced_molecules(
        entries,
        db_path=Path(db_path),
        vendor="COCONUT-2.0",
        source_id=src.id,
        source_url="https://coconut.naturalproducts.net/",
        license_spdx=src.spdx_license or src.license,
        in_stock=None,
        max_rows=max_rows,
    )


def parse_coconut_records(
    path: Path | str,
    *,
    dataset_id: str = "coconut_2_0",
    max_rows: int = 0,
) -> list[AssayfoldRecord]:
    """Convert COCONUT 2.0 natural-product rows into sealed `AssayfoldRecord(kind='molecule')` objects."""
    src = SOURCE_BY_ID["coconut"]
    lic = src.spdx_license or src.license
    out: list[AssayfoldRecord] = []
    for cid, mol_name, nm in _iter_coconut_rows(Path(path), max_rows=max_rows):
        mol = Molecule(
            smiles=nm.canonical_smiles,
            inchikey=nm.inchikey,
            name=mol_name or cid,
            formula=nm.formula,
        )
        prov = Provenance(
            source_url=f"https://coconut.naturalproducts.net/compounds/{cid}",
            used=[src.id],
        )
        rec = AssayfoldRecord(
            record_id=make_assayfold_record_id(src.id, dataset_id, cid),
            source_id=src.id,
            dataset_id=dataset_id,
            native_id=cid,
            kind="molecule",
            payload=mol,
            license=lic,
            origin="real",
            provenance=prov,
            tags=["corpus:coconut", "natural_product", f"license:{lic}"],
        ).sealed()
        out.append(rec)
    return out


# ---------- 5. ORDerly Reactions (Parquet / CSV / TSV -> reactions.sqlite) ----------


def _compose_orderly_rxn(row_rxn: Any, row_reactants: Any, row_products: Any, row_agents: Any = None) -> Optional[str]:
    if row_rxn and str(row_rxn).strip():
        return str(row_rxn).strip()
    if row_reactants and row_products:
        lhs = str(row_reactants).strip()
        rhs = str(row_products).strip()
        if lhs and rhs:
            if row_agents and str(row_agents).strip():
                return f"{lhs}>{str(row_agents).strip()}>>{rhs}"
            return f"{lhs}>>{rhs}"
    return None


def index_orderly_reactions(
    path: Path | str,
    db_path: Path | str = DEFAULT_REACTION_DB,
    *,
    max_rows: int = 0,
) -> int:
    """Index ORDerly `.parquet`, `.csv`, or `.tsv` reaction splits into `reactions.sqlite` (`CC-BY-4.0`)."""
    p = Path(path)

    def _pairs() -> Iterator[tuple[Optional[str], Optional[str]]]:
        if p.name.lower().endswith(".parquet"):
            cols = _read_parquet_columns(p)
            rxn_col = _pick_col(cols.keys(), ("rxn_str", "reaction_smiles", "rxn_smiles", "smiles"))
            react_col = _pick_col(cols.keys(), ("reactants", "reactant_smiles"))
            prod_col = _pick_col(cols.keys(), ("products", "product_smiles", "product"))
            agent_col = _pick_col(cols.keys(), ("solvents", "agents", "reagents"))
            id_col = _pick_col(cols.keys(), ("reaction_id", "ord_id", "id"))
            n = len(next(iter(cols.values()), []))
            for i in range(n):
                rxn = _compose_orderly_rxn(
                    cols[rxn_col][i] if rxn_col else None,
                    cols[react_col][i] if react_col else None,
                    cols[prod_col][i] if prod_col else None,
                    cols[agent_col][i] if agent_col else None,
                )
                ref = str(cols[id_col][i]) if (id_col and cols[id_col][i] is not None) else f"orderly:{i}"
                yield rxn, ref
            return

        for name, stream in _open_text_streams(p, suffixes=(".csv", ".tsv", ".txt")):
            delim = "\t" if name.lower().endswith(".tsv") else ","
            reader = csv.DictReader(stream, delimiter=delim)
            if not reader.fieldnames:
                continue
            rxn_col = _pick_col(reader.fieldnames, ("rxn_str", "reaction_smiles", "rxn_smiles", "smiles"))
            react_col = _pick_col(reader.fieldnames, ("reactants", "reactant_smiles"))
            prod_col = _pick_col(reader.fieldnames, ("products", "product_smiles", "product"))
            agent_col = _pick_col(reader.fieldnames, ("solvents", "agents", "reagents"))
            id_col = _pick_col(reader.fieldnames, ("reaction_id", "ord_id", "id"))
            for i, row in enumerate(reader):
                rxn = _compose_orderly_rxn(
                    row.get(rxn_col) if rxn_col else None,
                    row.get(react_col) if react_col else None,
                    row.get(prod_col) if prod_col else None,
                    row.get(agent_col) if agent_col else None,
                )
                ref = (row.get(id_col) if id_col else None) or f"orderly:{i}"
                yield rxn, ref

    return _index_reaction_pairs(
        _pairs(),
        db_path=Path(db_path),
        step_source="orderly_lookup",
        default_prefix="orderly",
        max_rows=max_rows,
    )


# ---------- 6. Syntheseus USPTO-50K & Retro* Inventory (MIT) ----------


def index_uspto_50k_reactions(
    path: Path | str,
    db_path: Path | str = DEFAULT_REACTION_DB,
    *,
    max_rows: int = 0,
) -> int:
    """Index Syntheseus USPTO-50K `.zip`, `.csv`, `.tsv`, or `.jsonl` reactions into `reactions.sqlite`."""
    p = Path(path)

    def _pairs() -> Iterator[tuple[Optional[str], Optional[str]]]:
        for name, stream in _open_text_streams(p, suffixes=(".csv", ".tsv", ".jsonl", ".smi", ".rsmi")):
            lname = name.lower()
            if lname.endswith(".jsonl"):
                for i, raw_line in enumerate(stream):
                    line = raw_line.strip()
                    if not line:
                        continue
                    obj = json.loads(line)
                    rxn = (
                        obj.get("rxn_smiles")
                        or obj.get("reactions")
                        or obj.get("reaction_smiles")
                        or _compose_orderly_rxn(None, obj.get("reactants"), obj.get("products"))
                    )
                    ref = obj.get("id") or f"uspto50k:{i}"
                    yield rxn, ref
            elif lname.endswith((".smi", ".rsmi")):
                for i, raw_line in enumerate(stream):
                    line = raw_line.strip()
                    if not line or line.startswith("#") or "ReactionSmiles" in line:
                        continue
                    parts = line.split()
                    yield parts[0], (parts[1] if len(parts) > 1 else f"uspto50k:{i}")
            else:
                delim = "\t" if lname.endswith(".tsv") else ","
                reader = csv.DictReader(stream, delimiter=delim)
                if not reader.fieldnames:
                    continue
                rxn_col = _pick_col(
                    reader.fieldnames,
                    ("reactions", "rxn_smiles", "reaction_smiles", "reactants>reagents>production", "smiles"),
                )
                react_col = _pick_col(reader.fieldnames, ("reactants", "reactant_smiles"))
                prod_col = _pick_col(reader.fieldnames, ("products", "prod_smiles", "product_smiles"))
                id_col = _pick_col(reader.fieldnames, ("id", "reaction_id", "patent_id"))
                for i, row in enumerate(reader):
                    rxn = _compose_orderly_rxn(
                        row.get(rxn_col) if rxn_col else None,
                        row.get(react_col) if react_col else None,
                        row.get(prod_col) if prod_col else None,
                    )
                    ref = (row.get(id_col) if id_col else None) or f"uspto50k:{i}"
                    yield rxn, ref

    return _index_reaction_pairs(
        _pairs(),
        db_path=Path(db_path),
        step_source="syntheseus_uspto_50k_lookup",
        default_prefix="uspto50k",
        max_rows=max_rows,
    )


def import_syntheseus_inventory(
    path: Path | str,
    db_path: Path | str = DEFAULT_STOCK_DB,
    *,
    max_rows: int = 0,
) -> int:
    """Import Syntheseus Retro* `ORIGINAL-inventory.csv[.gz]` into `stock.sqlite` (`MIT`)."""
    src = SOURCE_BY_ID["syntheseus_uspto_50k"]
    p = Path(path)

    def _iter() -> Iterator[tuple[NormalizedMol, Optional[str]]]:
        for name, stream in _open_text_streams(p, suffixes=(".csv", ".tsv", ".smi", ".txt")):
            delim = "\t" if name.lower().endswith(".tsv") else ","
            reader = csv.DictReader(stream, delimiter=delim)
            if not reader.fieldnames:
                continue
            smi_col = _pick_col(reader.fieldnames, ("smiles", "mol", "canonical_smiles"))
            id_col = _pick_col(reader.fieldnames, ("id", "catalog_id", "inchi_key", "inchikey"))
            if smi_col is None:
                continue
            for idx, row in enumerate(reader):
                raw_smi = row.get(smi_col)
                if not raw_smi:
                    continue
                nm = normalize_smiles(raw_smi, desalt=False)
                if nm is None:
                    continue
                cat_id = (row.get(id_col) or f"syntheseus-inv:{idx + 1}").strip() if id_col else f"syntheseus-inv:{idx + 1}"
                yield nm, cat_id

    return _insert_unpriced_molecules(
        _iter(),
        db_path=Path(db_path),
        vendor="Syntheseus-RetroStar",
        source_id=src.id,
        source_url="https://figshare.com/articles/dataset/25376728",
        license_spdx=src.spdx_license or src.license,
        in_stock=True,
        max_rows=max_rows,
    )


# ---------- 7. SynPlanner Open Data (building blocks + benchmark targets + USPTO reactions) ----------


def import_synplanner_building_blocks(
    path: Path | str,
    db_path: Path | str = DEFAULT_STOCK_DB,
    *,
    min_heavy_atoms: int = 2,
    max_rows: int = 0,
) -> int:
    """Import SynPlanner `building_blocks_*.smi` into `stock.sqlite`, skipping bare monoatomic ions."""
    src = SOURCE_BY_ID["synplanner_data"]
    p = Path(path)

    def _iter() -> Iterator[tuple[NormalizedMol, Optional[str]]]:
        idx = 0
        for _, stream in _open_text_streams(p, suffixes=(".smi", ".txt", ".csv", ".tsv")):
            for raw_line in stream:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if idx == 0 and parts[0].lower() in {"smiles", "mol"}:
                    continue
                nm = normalize_smiles(parts[0], desalt=False)
                if nm is None or nm.num_heavy_atoms < min_heavy_atoms:
                    continue
                idx += 1
                cat_id = parts[1] if len(parts) > 1 else f"synplanner-bb:{idx}"
                yield nm, cat_id

    return _insert_unpriced_molecules(
        _iter(),
        db_path=Path(db_path),
        vendor="SynPlanner-BB",
        source_id=src.id,
        source_url="https://huggingface.co/Laboratoire-De-Chemoinformatique/SynPlanner",
        license_spdx=src.spdx_license or src.license,
        in_stock=True,
        max_rows=max_rows,
    )


def load_synplanner_targets(path: Path | str, *, max_rows: int = 0) -> list[str]:
    """Load canonicalized benchmark target SMILES from SynPlanner `targets_with_sascore_*.smi`."""
    out: list[str] = []
    for _, stream in _open_text_streams(Path(path), suffixes=(".smi", ".txt")):
        for raw_line in stream:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            nm = normalize_smiles(parts[0], desalt=False)
            if nm is None:
                continue
            out.append(nm.canonical_smiles)
            if max_rows and len(out) >= max_rows:
                return out
    return out


def index_synplanner_reactions(
    path: Path | str,
    db_path: Path | str = DEFAULT_REACTION_DB,
    *,
    max_rows: int = 0,
) -> int:
    """Index SynPlanner standardized USPTO `.smi`, `.rsmi`, `.csv`, or `.zip` into `reactions.sqlite`."""
    p = Path(path)

    def _pairs() -> Iterator[tuple[Optional[str], Optional[str]]]:
        for name, stream in _open_text_streams(p, suffixes=(".smi", ".rsmi", ".csv", ".tsv", ".txt")):
            lname = name.lower()
            if lname.endswith((".csv", ".tsv")):
                delim = "\t" if lname.endswith(".tsv") else ","
                reader = csv.DictReader(stream, delimiter=delim)
                if not reader.fieldnames:
                    continue
                rxn_col = _pick_col(
                    reader.fieldnames,
                    ("reaction_smiles", "rxn_smiles", "reactions", "rsmi", "smiles"),
                )
                id_col = _pick_col(reader.fieldnames, ("id", "reaction_id", "patent_id"))
                if rxn_col is None:
                    continue
                for i, row in enumerate(reader):
                    yield row.get(rxn_col), (row.get(id_col) if id_col else f"synplanner:{i}")
            else:
                for i, raw_line in enumerate(stream):
                    line = raw_line.strip()
                    if not line or line.startswith("#") or ">>" not in line:
                        continue
                    parts = line.split()
                    yield parts[0], (parts[1] if len(parts) > 1 else f"synplanner:{i}")

    return _index_reaction_pairs(
        _pairs(),
        db_path=Path(db_path),
        step_source="synplanner_lookup",
        default_prefix="synplanner",
        max_rows=max_rows,
    )


def ingest_downloaded_artifact(
    artifact_key: str,
    file_path: Path | str,
    *,
    stock_db: Path | str = DEFAULT_STOCK_DB,
    reference_db: Path | str = DEFAULT_REFERENCE_COMPOUNDS_DB,
    reaction_db: Path | str = DEFAULT_REACTION_DB,
    max_rows: int = 0,
) -> int:
    """Dispatch a downloaded open-corpus file to its matching SQLite importer/indexer."""
    p = Path(file_path)
    if artifact_key == "paroutes_stock_n1":
        return import_paroutes_stock(p, db_path=stock_db, subset="n1", max_rows=max_rows)
    if artifact_key == "paroutes_stock_n5":
        return import_paroutes_stock(p, db_path=stock_db, subset="n5", max_rows=max_rows)
    if artifact_key in {"paroutes_routes_n1", "paroutes_routes_n5"}:
        return index_paroutes_routes(p, db_path=reaction_db, max_rows=max_rows)
    if artifact_key == "paroutes_reactions":
        return index_paroutes_reactions_csv(p, db_path=reaction_db, max_rows=max_rows)
    if artifact_key == "chembl_37_chemreps":
        return import_chembl_chemreps(p, db_path=reference_db, release="chembl_37", max_rows=max_rows)
    if artifact_key == "surechembl_compounds_parquet":
        return import_surechembl_compounds(p, db_path=reference_db, max_rows=max_rows)
    if artifact_key in {"coconut_csv_lite_2026_11", "coconut_2_0_zenodo_csv"}:
        return import_coconut_compounds(p, db_path=reference_db, max_rows=max_rows)
    if artifact_key in {"orderly_retro_test", "orderly_retro_train"}:
        return index_orderly_reactions(p, db_path=reaction_db, max_rows=max_rows)
    if artifact_key == "syntheseus_uspto_50k":
        return index_uspto_50k_reactions(p, db_path=reaction_db, max_rows=max_rows)
    if artifact_key == "syntheseus_retro_star_inventory":
        return import_syntheseus_inventory(p, db_path=stock_db, max_rows=max_rows)
    if artifact_key == "synplanner_building_blocks":
        return import_synplanner_building_blocks(p, db_path=stock_db, max_rows=max_rows)
    if artifact_key == "synplanner_targets_sa_1_5_2_5":
        return len(load_synplanner_targets(p, max_rows=max_rows))
    if artifact_key == "synplanner_uspto_standardized":
        return index_synplanner_reactions(p, db_path=reaction_db, max_rows=max_rows)
    raise KeyError(f"Unsupported artifact_key {artifact_key!r}")
