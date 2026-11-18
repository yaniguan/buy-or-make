"""Stock and vendor offer backends (local files/SQLite + live PubChem PUG-REST/PUG-View)."""

from __future__ import annotations

import csv
import gzip
import json
import sqlite3
import time
from pathlib import Path
from typing import Iterable, Optional

import httpx

from .cache import ResponseCache
from .canonical import normalize_smiles
from .schema import Offer, Quantity, utc_now

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest"


def _parse_pack_size(val: Optional[str], unit: Optional[str]) -> Optional[Quantity]:
    if val is None or str(val).strip() == "":
        return None
    try:
        num = float(val)
    except (TypeError, ValueError):
        return None
    u = (unit or "1").strip() or "1"
    return Quantity(value=num, unit=u, raw_unit_str=f"{val} {u}".strip())


def _parse_optional_float(val) -> Optional[float]:
    if val is None or str(val).strip() == "":
        return None
    return float(val)


def _parse_optional_int(val) -> Optional[int]:
    if val is None or str(val).strip() == "":
        return None
    return int(val)


def _parse_optional_bool(val) -> Optional[bool]:
    if val is None or str(val).strip() == "":
        return None
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower()
    if s in {"1", "true", "yes", "y"}:
        return True
    if s in {"0", "false", "no", "n"}:
        return False
    return None


class LocalStockBackend:
    """Lookup purchasable compounds and optional catalog offers from local files or SQLite."""

    def __init__(self, paths: Optional[Iterable[Path | str]] = None) -> None:
        self._offers_by_ik: dict[str, list[Offer]] = {}
        self._sqlite_conns: list[sqlite3.Connection] = []
        if paths:
            for p in paths:
                self.load_path(Path(p))

    def load_path(self, path: Path) -> None:
        if not path.exists():
            return
        name = path.name.lower()
        if name.endswith((".sqlite", ".db")):
            self._sqlite_conns.append(sqlite3.connect(str(path)))
        elif name.endswith(".jsonl"):
            self._load_jsonl(path)
        elif name.endswith((".tsv", ".csv")):
            self._load_delimited(path, delimiter="\t" if name.endswith(".tsv") else ",")
        elif name.endswith((".smi", ".smi.gz", ".txt")):
            self._load_smi(path)
        elif name.endswith((".hdf5", ".h5")):
            self._load_hdf5_stock(path)

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
                pack = _parse_pack_size(row.get("pack_value") or row.get("pack_size"), row.get("pack_unit"))
                offer = Offer(
                    inchikey=ik,
                    vendor=(row.get("vendor") or "local_stock").strip(),
                    catalog_id=(row.get("catalog_id") or "").strip() or None,
                    pack_size=pack,
                    purity=(row.get("purity") or "").strip() or None,
                    price=price,
                    currency=currency,
                    price_basis="listed_catalog" if price is not None else "unknown",
                    lead_time_days=_parse_optional_int(row.get("lead_time_days")),
                    in_stock=_parse_optional_bool(row.get("in_stock")),
                    source_id=(row.get("source_id") or "local_stock").strip(),
                    source_url=(row.get("source_url") or path.as_uri()).strip(),
                    fetched_utc=(row.get("fetched_utc") or utc_now()).strip(),
                    license=(row.get("license") or "custom/terms").strip(),
                )
                self._add_offer(offer)

    def _load_smi(self, path: Path) -> None:
        opener = gzip.open if path.name.lower().endswith(".gz") else open
        ts = utc_now()
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
                        vendor=path.stem,
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
        """Load AiZynthFinder-style HDF5 stock if PyTables/pandas or h5py is installed."""
        ts = utc_now()
        try:
            import pandas as pd  # type: ignore

            df = pd.read_hdf(path, key="table")
            col = "inchi_key" if "inchi_key" in df.columns else "inchikey"
            for ik in df[col].dropna().astype(str):
                self._add_offer(
                    Offer(
                        inchikey=ik.strip(),
                        vendor="zinc_stock_hdf5",
                        price=None,
                        currency=None,
                        price_basis="unknown",
                        source_id="aizynthfinder_zinc_stock",
                        source_url=path.resolve().as_uri(),
                        fetched_utc=ts,
                        license="custom/terms",
                    )
                )
        except Exception:
            return

    def lookup(self, inchikey: str) -> list[Offer]:
        offers = list(self._offers_by_ik.get(inchikey, []))
        for conn in self._sqlite_conns:
            rows = conn.execute(
                "SELECT vendor, catalog_id, source_id, source_url, fetched_utc, license "
                "FROM stock WHERE inchikey = ?",
                (inchikey,),
            ).fetchall()
            for vendor, cat_id, src_id, src_url, fetched_utc, lic in rows:
                offers.append(
                    Offer(
                        inchikey=inchikey,
                        vendor=vendor,
                        catalog_id=cat_id,
                        price=None,
                        currency=None,
                        price_basis="unknown",
                        source_id=src_id,
                        source_url=src_url,
                        fetched_utc=fetched_utc,
                        license=lic,
                    )
                )
        return offers

    def close(self) -> None:
        for conn in self._sqlite_conns:
            conn.close()


class PubChemVendorBackend:
    """Query PubChem PUG-REST + PUG-View Chemical Vendors with rate limiting and caching."""

    def __init__(
        self,
        cache: ResponseCache,
        *,
        offline: bool = False,
        min_interval_sec: float = 0.22,
        max_offers: int = 25,
        client: Optional[httpx.Client] = None,
    ) -> None:
        self.cache = cache
        self.offline = offline
        self.min_interval_sec = min_interval_sec
        self.max_offers = max_offers
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=15.0,
            headers={"User-Agent": "buy-or-make/0.1 (open chemistry availability check)"},
        )
        self._last_req_ts = 0.0

    def _throttle(self) -> None:
        if self.min_interval_sec <= 0:
            return
        elapsed = time.monotonic() - self._last_req_ts
        if elapsed < self.min_interval_sec:
            time.sleep(self.min_interval_sec - elapsed)
        self._last_req_ts = time.monotonic()

    def _get_json(self, url: str) -> Optional[dict]:
        self._throttle()
        resp = self._client.get(url)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        return resp.json()

    def lookup(self, inchikey: str, canonical_smiles: Optional[str] = None) -> list[Offer]:
        cache_key = f"pubchem:vendors:v1:{inchikey}"
        cached = self.cache.get(cache_key)
        if cached is not None:
            payload, fetched_utc = cached
            return self._payload_to_offers(inchikey, payload, fetched_utc)

        if self.offline:
            return []

        cid = self._fetch_cid(inchikey, canonical_smiles)
        if cid is None:
            payload = {"cid": None, "vendors": []}
            self.cache.set(cache_key, "pubchem_pug", payload)
            return []

        vendors = self._fetch_chemical_vendors(cid)
        payload = {"cid": cid, "vendors": vendors}
        fetched_utc = self.cache.set(cache_key, "pubchem_pug", payload)
        return self._payload_to_offers(inchikey, payload, fetched_utc)

    def _fetch_cid(self, inchikey: str, canonical_smiles: Optional[str]) -> Optional[int]:
        try:
            data = self._get_json(f"{PUBCHEM_BASE}/pug/compound/inchikey/{inchikey}/cids/JSON")
        except httpx.HTTPError:
            return None
        if data and "IdentifierList" in data and data["IdentifierList"].get("CID"):
            cids = [int(c) for c in data["IdentifierList"]["CID"] if int(c) > 0]
            if cids:
                return cids[0]
        if canonical_smiles:
            try:
                resp = self._client.post(
                    f"{PUBCHEM_BASE}/pug/compound/smiles/cids/JSON",
                    data={"smiles": canonical_smiles},
                )
                self._last_req_ts = time.monotonic()
                if resp.status_code == 200:
                    d2 = resp.json()
                    cids = [int(c) for c in d2.get("IdentifierList", {}).get("CID", []) if int(c) > 0]
                    if cids:
                        return cids[0]
            except httpx.HTTPError:
                return None
        return None

    def _fetch_chemical_vendors(self, cid: int) -> list[dict]:
        try:
            data = self._get_json(f"{PUBCHEM_BASE}/pug_view/categories/compound/{cid}/JSON")
        except httpx.HTTPError:
            return []
        if not data:
            return []
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
                    return vendors
        return vendors

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
