"""Stock and vendor offer backends (local files/SQLite/HDF5 + live PubChem PUG-REST/PUG-View)."""

from __future__ import annotations

import csv
import gzip
import json
import sqlite3
import time
import warnings
from pathlib import Path
from typing import Any, Iterable, Optional

import httpx
import numpy as np

from .cache import ResponseCache
from .canonical import normalize_smiles
from .schema import Offer, Quantity, utc_now

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest"


def _parse_pack_size(val: Optional[Any], unit: Optional[str]) -> Optional[Quantity]:
    if val is None or str(val).strip() == "":
        return None
    try:
        num = float(val)
    except (TypeError, ValueError):
        return None
    if num <= 0:
        return None
    u = (unit or "1").strip() or "1"
    return Quantity(value=num, unit=u, raw_unit_str=f"{val} {u}".strip())


def _parse_optional_float(val: Any) -> Optional[float]:
    if val is None or isinstance(val, bool) or str(val).strip() == "":
        return None
    return float(val)


def _parse_optional_int(val: Any) -> Optional[int]:
    if val is None or isinstance(val, bool) or str(val).strip() == "":
        return None
    return int(val)


def _parse_optional_bool(val: Any) -> Optional[bool]:
    if val is None or str(val).strip() == "":
        return None
    if isinstance(val, bool):
        return val
    if isinstance(val, int):
        return bool(val)
    s = str(val).strip().lower()
    if s in {"1", "true", "yes", "y"}:
        return True
    if s in {"0", "false", "no", "n"}:
        return False
    return None


def _decode_h5_str(val: Any) -> str:
    if isinstance(val, (bytes, bytearray, np.bytes_)):
        return bytes(val).decode("utf-8", errors="replace").strip()
    return str(val).strip()


def _extract_inchikeys_via_h5py(h5_obj: Any) -> Optional[list[str]]:
    """Recursively locate an `inchi_key` / `inchikey` column or dataset in an open h5py File/Group."""
    # PyTables Table compound dataset (e.g. /table/table with structured dtype)
    if hasattr(h5_obj, "dtype") and getattr(h5_obj.dtype, "names", None):
        names = h5_obj.dtype.names
        for key in ("inchi_key", "inchikey"):
            if key in names:
                col_vals = h5_obj[key]
                return [s for raw in col_vals if (s := _decode_h5_str(raw))]
        return None

    # Group / File mapping containing sub-groups or 1D datasets
    if hasattr(h5_obj, "keys"):
        keys = set(h5_obj.keys())
        for key in ("inchi_key", "inchikey"):
            if key in keys and hasattr(h5_obj[key], "shape"):
                ds = h5_obj[key][:]
                return [s for raw in ds if (s := _decode_h5_str(raw))]
        for sub_name in h5_obj.keys():
            found = _extract_inchikeys_via_h5py(h5_obj[sub_name])
            if found is not None:
                return found
    return None


class LocalStockBackend:
    """Lookup purchasable compounds and optional catalog offers from local files or SQLite."""

    def __init__(self, paths: Optional[Iterable[Path | str]] = None) -> None:
        self._offers_by_ik: dict[str, list[Offer]] = {}
        self._sqlite_tables: list[tuple[sqlite3.Connection, Path, tuple[str, ...]]] = []
        if paths:
            for p in paths:
                self.load_path(Path(p))

    def load_path(self, path: Path) -> None:
        if not path.exists():
            return
        name = path.name.lower()
        if name.endswith((".sqlite", ".db")):
            self._open_sqlite(path)
        elif name.endswith(".jsonl"):
            self._load_jsonl(path)
        elif name.endswith((".tsv", ".csv")):
            self._load_delimited(path, delimiter="\t" if name.endswith(".tsv") else ",")
        elif name.endswith((".smi", ".smi.gz", ".txt")):
            self._load_smi(path)
        elif name.endswith((".hdf5", ".h5")):
            self._load_hdf5_stock(path)

    def _open_sqlite(self, path: Path) -> None:
        conn = sqlite3.connect(str(path))
        cur = conn.execute("PRAGMA table_info(stock)")
        cols = tuple(str(row[1]) for row in cur.fetchall())
        if not cols:
            conn.close()
            warnings.warn(f"SQLite stock file '{path}' has no 'stock' table; skipping.", stacklevel=2)
            return
        self._sqlite_tables.append((conn, path, cols))

    def _add_offer(self, offer: Offer) -> None:
        self._offers_by_ik.setdefault(offer.inchikey, []).append(offer)

    def _load_jsonl(self, path: Path) -> None:
        with path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                row = json.loads(line)
                ik = row.get("inchikey")
                if not ik and row.get("smiles"):
                    nm = normalize_smiles(row["smiles"])
                    ik = nm.inchikey if nm else None
                if not ik:
                    continue
                row["inchikey"] = ik
                self._add_offer(Offer.model_validate(row))

    def _load_delimited(self, path: Path, delimiter: str) -> None:
        with path.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                ik = (row.get("inchikey") or "").strip()
                if not ik and row.get("smiles"):
                    nm = normalize_smiles(row["smiles"])
                    ik = nm.inchikey if nm else ""
                if not ik:
                    continue
                price = _parse_optional_float(row.get("price"))
                currency = (row.get("currency") or "").strip() or None
                if price is None:
                    currency = None
                elif currency is None:
                    currency = "USD"
                pack = _parse_pack_size(
                    row.get("pack_value") or row.get("pack_size"), row.get("pack_unit")
                )
                basis_raw = (row.get("price_basis") or "").strip()
                price_basis = (
                    basis_raw
                    if basis_raw in {"listed_catalog", "unknown"} and not (basis_raw == "listed_catalog" and price is None)
                    else ("listed_catalog" if price is not None else "unknown")
                )
                offer = Offer(
                    inchikey=ik,
                    vendor=(row.get("vendor") or "local_stock").strip(),
                    catalog_id=(row.get("catalog_id") or "").strip() or None,
                    pack_size=pack,
                    purity=(row.get("purity") or "").strip() or None,
                    price=price,
                    currency=currency,
                    price_basis=price_basis,  # type: ignore[arg-type]
                    lead_time_days=_parse_optional_int(row.get("lead_time_days")),
                    in_stock=_parse_optional_bool(row.get("in_stock")),
                    source_id=(row.get("source_id") or "local_stock").strip(),
                    source_url=(row.get("source_url") or path.resolve().as_uri()).strip(),
                    fetched_utc=(row.get("fetched_utc") or utc_now()).strip(),
                    license=(row.get("license") or "custom/terms").strip(),
                )
                self._add_offer(offer)

    def _load_smi(self, path: Path) -> None:
        opener = gzip.open if path.name.lower().endswith(".gz") else open
        ts = utc_now()
        vendor_name = path.name[:-7] if path.name.lower().endswith(".smi.gz") else path.stem
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split()
                if not parts or parts[0].startswith("#"):
                    continue
                nm = normalize_smiles(parts[0])
                if nm is None:
                    continue
                cat_id = parts[1] if len(parts) > 1 else None
                self._add_offer(
                    Offer(
                        inchikey=nm.inchikey,
                        vendor=vendor_name,
                        catalog_id=cat_id,
                        price=None,
                        currency=None,
                        price_basis="unknown",
                        source_id="local_stock",
                        source_url=path.resolve().as_uri(),
                        fetched_utc=ts,
                        license="custom/terms",
                    )
                )

    def _load_hdf5_stock(self, path: Path) -> None:
        """Load AiZynthFinder-style or generic HDF5 stock via pandas/PyTables or h5py fallback."""
        ts = utc_now()
        errors: list[str] = []
        inchikeys: Optional[list[str]] = None

        try:
            import pandas as pd  # type: ignore

            df = pd.read_hdf(path, key="table")
            col = "inchi_key" if "inchi_key" in df.columns else ("inchikey" if "inchikey" in df.columns else None)
            if col is None:
                errors.append(f"pandas.read_hdf found columns {list(df.columns)}, expected 'inchi_key' or 'inchikey'")
            else:
                inchikeys = [s.strip() for s in df[col].dropna().astype(str) if s.strip()]
        except Exception as exc:
            errors.append(f"pandas/tables: {type(exc).__name__}: {exc}")

        if inchikeys is None:
            try:
                import h5py  # type: ignore

                with h5py.File(path, "r") as hf:
                    inchikeys = _extract_inchikeys_via_h5py(hf)
                    if inchikeys is None:
                        errors.append("h5py: no 'inchi_key' or 'inchikey' dataset found in HDF5 hierarchy")
            except Exception as exc:
                errors.append(f"h5py: {type(exc).__name__}: {exc}")

        if inchikeys is None:
            detail = "; ".join(errors)
            msg = (
                f"Failed to read HDF5 stock file '{path}'. Install HDF5 support via "
                f"`pip install 'buy-or-make[hdf5]'` (h5py or PyTables) and verify dataset layout. "
                f"Details: {detail}"
            )
            warnings.warn(msg, UserWarning, stacklevel=2)
            raise RuntimeError(msg)

        src_uri = path.resolve().as_uri()
        for ik in inchikeys:
            self._add_offer(
                Offer(
                    inchikey=ik,
                    vendor="zinc_stock_hdf5",
                    price=None,
                    currency=None,
                    price_basis="unknown",
                    source_id="aizynthfinder_zinc_stock",
                    source_url=src_uri,
                    fetched_utc=ts,
                    license="custom/terms",
                )
            )

    def lookup(self, inchikey: str) -> list[Offer]:
        offers = list(self._offers_by_ik.get(inchikey, []))
        desired_cols = (
            "vendor",
            "catalog_id",
            "source_id",
            "source_url",
            "fetched_utc",
            "license",
            "price",
            "currency",
            "pack_value",
            "pack_size",
            "pack_unit",
            "purity",
            "lead_time_days",
            "in_stock",
            "price_basis",
        )
        for conn, db_path, available_cols in self._sqlite_tables:
            select_cols = [c for c in desired_cols if c in available_cols]
            if not select_cols:
                continue
            col_sql = ", ".join(f'"{c}"' for c in select_cols)
            rows = conn.execute(
                f'SELECT {col_sql} FROM stock WHERE "inchikey" = ?',
                (inchikey,),
            ).fetchall()
            for raw_vals in rows:
                row = dict(zip(select_cols, raw_vals))
                price = _parse_optional_float(row.get("price"))
                currency = (str(row["currency"]).strip() if row.get("currency") is not None else "") or None
                if price is None:
                    currency = None
                elif currency is None:
                    currency = "USD"
                pack_val = row.get("pack_value") if row.get("pack_value") is not None else row.get("pack_size")
                pack_unit = str(row["pack_unit"]).strip() if row.get("pack_unit") is not None else None
                pack = _parse_pack_size(pack_val, pack_unit)
                basis_raw = str(row.get("price_basis") or "").strip()
                price_basis = (
                    basis_raw
                    if basis_raw in {"listed_catalog", "unknown"} and not (basis_raw == "listed_catalog" and price is None)
                    else ("listed_catalog" if price is not None else "unknown")
                )
                offers.append(
                    Offer(
                        inchikey=inchikey,
                        vendor=str(row.get("vendor") or "local_stock"),
                        catalog_id=(str(row["catalog_id"]).strip() if row.get("catalog_id") is not None else None) or None,
                        pack_size=pack,
                        purity=(str(row["purity"]).strip() if row.get("purity") is not None else None) or None,
                        price=price,
                        currency=currency,
                        price_basis=price_basis,  # type: ignore[arg-type]
                        lead_time_days=_parse_optional_int(row.get("lead_time_days")),
                        in_stock=_parse_optional_bool(row.get("in_stock")),
                        source_id=str(row.get("source_id") or "local_stock"),
                        source_url=str(row.get("source_url") or db_path.resolve().as_uri()),
                        fetched_utc=str(row.get("fetched_utc") or utc_now()),
                        license=str(row.get("license") or "custom/terms"),
                    )
                )
        return offers

    def close(self) -> None:
        for conn, _, _ in self._sqlite_tables:
            conn.close()
        self._sqlite_tables.clear()


class PubChemVendorBackend:
    """Query PubChem PUG-REST + PUG-View Chemical Vendors with rate limiting, retries, and safe caching."""

    def __init__(
        self,
        cache: ResponseCache,
        *,
        offline: bool = False,
        min_interval_sec: float = 0.22,
        max_offers: int = 25,
        max_retries: int = 2,
        backoff_base_sec: float = 0.20,
        client: Optional[httpx.Client] = None,
    ) -> None:
        self.cache = cache
        self.offline = offline
        self.min_interval_sec = min_interval_sec
        self.max_offers = max_offers
        self.max_retries = max(0, int(max_retries))
        self.backoff_base_sec = max(0.0, float(backoff_base_sec))
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=15.0,
            headers={"User-Agent": "buy-or-make/0.1 (open chemistry availability check)"},
        )
        self._last_req_ts = 0.0

    def _throttle(self) -> None:
        if self.min_interval_sec <= 0:
            self._last_req_ts = time.monotonic()
            return
        elapsed = time.monotonic() - self._last_req_ts
        if elapsed < self.min_interval_sec:
            time.sleep(self.min_interval_sec - elapsed)
        self._last_req_ts = time.monotonic()

    def _request_json(
        self,
        method: str,
        url: str,
        *,
        data: Optional[dict[str, str]] = None,
    ) -> tuple[Optional[dict], bool]:
        """Execute a throttled HTTP request with retry on transient errors.

        Returns `(json_data, definitive_not_found)`:
        - `(dict, False)`: HTTP 200 valid JSON response (safe to cache).
        - `(None, True)`: Definitive negative (`404` / `400` `PUGREST.NotFound`, safe to cache).
        - `(None, False)`: Transient transport/429/5xx failure after retries (MUST NOT cache!).
        """
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                if method == "GET":
                    resp = self._client.get(url)
                else:
                    resp = self._client.post(url, data=data)
                self._last_req_ts = time.monotonic()
            except httpx.TransportError:
                if attempt < self.max_retries and self.backoff_base_sec > 0:
                    time.sleep(self.backoff_base_sec * (2**attempt))
                continue
            except httpx.HTTPError:
                return None, False

            if resp.status_code in {400, 404}:
                return None, True
            if resp.status_code == 429 or resp.status_code >= 500:
                if attempt < self.max_retries and self.backoff_base_sec > 0:
                    time.sleep(self.backoff_base_sec * (2**attempt))
                continue
            if 200 <= resp.status_code < 300:
                try:
                    return resp.json(), False
                except ValueError:
                    return None, False
            # Other 4xx statuses (e.g. 401/403 proxy errors) are not cached as negative compounds
            return None, False
        return None, False

    def _get_json(self, url: str) -> Optional[dict]:
        data, _ = self._request_json("GET", url)
        return data

    def lookup(self, inchikey: str, canonical_smiles: Optional[str] = None) -> list[Offer]:
        cache_key = f"pubchem:vendors:v1:{inchikey}"
        cached = self.cache.get(cache_key)
        if cached is not None:
            payload, fetched_utc = cached
            return self._payload_to_offers(inchikey, payload, fetched_utc)

        if self.offline:
            return []

        cid, cid_cacheable = self._fetch_cid(inchikey, canonical_smiles)
        if cid is None:
            if cid_cacheable:
                payload = {"cid": None, "vendors": []}
                self.cache.set(cache_key, "pubchem_pug", payload)
            return []

        vendors, vendors_cacheable = self._fetch_chemical_vendors(cid)
        if not vendors_cacheable:
            return []
        payload = {"cid": cid, "vendors": vendors}
        fetched_utc = self.cache.set(cache_key, "pubchem_pug", payload)
        return self._payload_to_offers(inchikey, payload, fetched_utc)

    def _fetch_cid(
        self,
        inchikey: str,
        canonical_smiles: Optional[str],
    ) -> tuple[Optional[int], bool]:
        data, not_found = self._request_json(
            "GET",
            f"{PUBCHEM_BASE}/pug/compound/inchikey/{inchikey}/cids/JSON",
        )
        transient_failure = data is None and not not_found
        if data and "IdentifierList" in data and data["IdentifierList"].get("CID"):
            cids = [int(c) for c in data["IdentifierList"]["CID"] if int(c) > 0]
            if cids:
                return cids[0], True

        if canonical_smiles:
            d2, nf2 = self._request_json(
                "POST",
                f"{PUBCHEM_BASE}/pug/compound/smiles/cids/JSON",
                data={"smiles": canonical_smiles},
            )
            if d2 is None and not nf2:
                transient_failure = True
            elif d2 and "IdentifierList" in d2 and d2["IdentifierList"].get("CID"):
                cids = [int(c) for c in d2["IdentifierList"]["CID"] if int(c) > 0]
                if cids:
                    return cids[0], True

        return None, (not transient_failure)

    def _fetch_chemical_vendors(self, cid: int) -> tuple[list[dict], bool]:
        data, not_found = self._request_json(
            "GET",
            f"{PUBCHEM_BASE}/pug_view/categories/compound/{cid}/JSON",
        )
        if data is None:
            # 404 on PUG-View means compound has no PUG-View vendor categories -> cacheable empty list.
            # Transient error (not_found=False) -> NOT cacheable!
            return [], not_found
        cats = data.get("SourceCategories", {}).get("Categories", [])
        vendors: list[dict] = []
        seen: set[tuple[str, Optional[str]]] = set()
        for cat in cats:
            if cat.get("Category") != "Chemical Vendors":
                continue
            for src in cat.get("Sources", []):
                name = (src.get("SourceName") or "").strip()
                reg_id = (src.get("RegistryID") or "").strip() or None
                if not name or (name, reg_id) in seen:
                    continue
                seen.add((name, reg_id))
                vendors.append(
                    {
                        "vendor": name,
                        "catalog_id": reg_id,
                        "source_url": src.get("SourceRecordURL")
                        or src.get("SourceDetail")
                        or f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=Chemical-Vendors",
                    }
                )
                if len(vendors) >= self.max_offers:
                    return vendors, True
        return vendors, True

    def _payload_to_offers(self, inchikey: str, payload: dict, fetched_utc: str) -> list[Offer]:
        offers: list[Offer] = []
        for item in payload.get("vendors", []):
            offers.append(
                Offer(
                    inchikey=inchikey,
                    vendor=item["vendor"],
                    catalog_id=item.get("catalog_id"),
                    price=None,
                    currency=None,
                    price_basis="unknown",
                    in_stock=None,
                    source_id="pubchem_pug",
                    source_url=item["source_url"],
                    fetched_utc=fetched_utc,
                    license="custom/terms",
                )
            )
        return offers

    def close(self) -> None:
        if self._owns_client:
            self._client.close()
