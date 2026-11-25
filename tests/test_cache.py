"""Tests for PubChem caching, transient error guards, throttling, TTL, and offline mode."""

from __future__ import annotations

from pathlib import Path

import httpx

from buyormake import cost_aware_rank
from buyormake.cache import ResponseCache
from buyormake.canonical import normalize_smiles
from buyormake.stock import PubChemVendorBackend


def test_pubchem_cache_and_offline_mode(tmp_path: Path) -> None:
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    nm = normalize_smiles(aspirin)
    assert nm is not None

    calls = {"count": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        calls["count"] += 1
        url = str(request.url)
        if f"/inchikey/{nm.inchikey}/cids/JSON" in url:
            return httpx.Response(200, json={"IdentifierList": {"CID": [2244]}})
        if "/categories/compound/2244/JSON" in url:
            return httpx.Response(
                200,
                json={
                    "SourceCategories": {
                        "Categories": [
                            {
                                "Category": "Chemical Vendors",
                                "Sources": [
                                    {
                                        "SID": 123,
                                        "SourceName": "MockVendorChem",
                                        "RegistryID": "MV-2244",
                                        "SourceRecordURL": "https://mockvendor.example/MV-2244",
                                    }
                                ],
                            }
                        ]
                    }
                },
            )
        return httpx.Response(404, text="Not Found")

    cache_path = tmp_path / "cache.sqlite"
    mock_client = httpx.Client(transport=httpx.MockTransport(handler))

    # 1. Cold online call against mock transport populates cache
    with ResponseCache(cache_path) as cache:
        backend = PubChemVendorBackend(cache, offline=False, min_interval_sec=0.0, client=mock_client)
        offers = backend.lookup(nm.inchikey, nm.canonical_smiles)
        assert len(offers) == 1
        assert offers[0].vendor == "MockVendorChem"
        assert offers[0].catalog_id == "MV-2244"
        assert offers[0].price is None
        assert offers[0].price_basis == "unknown"
        assert calls["count"] == 2

    # 2. Offline call must hit SQLite cache and NOT touch HTTP transport
    def fail_handler(request: httpx.Request) -> httpx.Response:
        raise AssertionError(f"Network call attempted in offline mode: {request.url}")

    fail_client = httpx.Client(transport=httpx.MockTransport(fail_handler))
    results = cost_aware_rank(
        [aspirin, "CCCCCC(=O)NCCCC"],  # aspirin cached, second uncached
        max_steps=1,
        backends=("pubchem",),
        offline=True,
        cache_path=cache_path,
        stock_paths=[],
        reaction_db_path=tmp_path / "missing_rxn.sqlite",
        scscore_weights=tmp_path / "missing_sc.json.gz",
        http_client=fail_client,
    )
    assert results[0].status == "buyable"
    assert results[0].best_price is None
    assert results[0].evidence == ["https://mockvendor.example/MV-2244"]
    assert results[1].status == "unknown"


def test_transient_503_and_timeout_never_poison_cache(tmp_path: Path) -> None:
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    nm = normalize_smiles(aspirin)
    assert nm is not None
    cache_key = f"pubchem:vendors:v1:{nm.inchikey}"

    # Case A: 503 Server Busy on CID lookup -> retries exhausted -> NOT cached
    def busy_cid_handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(503, text="PUGREST.ServerBusy")

    with ResponseCache(tmp_path / "c1.sqlite") as cache:
        backend = PubChemVendorBackend(
            cache,
            min_interval_sec=0.0,
            max_retries=1,
            backoff_base_sec=0.0,
            client=httpx.Client(transport=httpx.MockTransport(busy_cid_handler)),
        )
        assert backend.lookup(nm.inchikey, nm.canonical_smiles) == []
        assert cache.get(cache_key) is None

    # Case B: ConnectTimeout on CID lookup -> NOT cached
    def timeout_handler(request: httpx.Request) -> httpx.Response:
        raise httpx.ConnectTimeout("connection timed out", request=request)

    with ResponseCache(tmp_path / "c2.sqlite") as cache:
        backend = PubChemVendorBackend(
            cache,
            min_interval_sec=0.0,
            max_retries=1,
            backoff_base_sec=0.0,
            client=httpx.Client(transport=httpx.MockTransport(timeout_handler)),
        )
        assert backend.lookup(nm.inchikey, nm.canonical_smiles) == []
        assert cache.get(cache_key) is None

    # Case C: CID succeeds (200 -> 2244), but PUG-View returns 503 -> aspirin MUST NOT be cached as empty!
    def pugview_503_handler(request: httpx.Request) -> httpx.Response:
        if "/inchikey/" in str(request.url):
            return httpx.Response(200, json={"IdentifierList": {"CID": [2244]}})
        return httpx.Response(503, text="Service Unavailable")

    with ResponseCache(tmp_path / "c3.sqlite") as cache:
        backend = PubChemVendorBackend(
            cache,
            min_interval_sec=0.0,
            max_retries=1,
            backoff_base_sec=0.0,
            client=httpx.Client(transport=httpx.MockTransport(pugview_503_handler)),
        )
        assert backend.lookup(nm.inchikey, nm.canonical_smiles) == []
        assert cache.get(cache_key) is None

    # Case D: Transient 429 on first PUG-View attempt recovers on retry -> cached with vendor!
    attempts = {"pugview": 0}

    def retry_recovers_handler(request: httpx.Request) -> httpx.Response:
        if "/inchikey/" in str(request.url):
            return httpx.Response(200, json={"IdentifierList": {"CID": [2244]}})
        attempts["pugview"] += 1
        if attempts["pugview"] == 1:
            return httpx.Response(429, text="Too Many Requests")
        return httpx.Response(
            200,
            json={
                "SourceCategories": {
                    "Categories": [
                        {
                            "Category": "Chemical Vendors",
                            "Sources": [{"SourceName": "RecoveredVendor", "RegistryID": "RV-1"}],
                        }
                    ]
                }
            },
        )

    with ResponseCache(tmp_path / "c4.sqlite") as cache:
        backend = PubChemVendorBackend(
            cache,
            min_interval_sec=0.0,
            max_retries=2,
            backoff_base_sec=0.0,
            client=httpx.Client(transport=httpx.MockTransport(retry_recovers_handler)),
        )
        offers = backend.lookup(nm.inchikey, nm.canonical_smiles)
        assert len(offers) == 1
        assert offers[0].vendor == "RecoveredVendor"
        assert cache.get(cache_key) is not None

    # Case E: Definitive 404 on both InChIKey GET and SMILES POST IS cached as negative result
    def not_found_handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(404, json={"Fault": {"Code": "PUGREST.NotFound"}})

    with ResponseCache(tmp_path / "c5.sqlite") as cache:
        backend = PubChemVendorBackend(
            cache,
            min_interval_sec=0.0,
            max_retries=1,
            backoff_base_sec=0.0,
            client=httpx.Client(transport=httpx.MockTransport(not_found_handler)),
        )
        assert backend.lookup(nm.inchikey, nm.canonical_smiles) == []
        cached_neg = cache.get(cache_key)
        assert cached_neg is not None
        assert cached_neg[0] == {"cid": None, "vendors": []}


def test_smiles_post_fallback_and_throttle_on_every_method(tmp_path: Path, monkeypatch) -> None:
    nm = normalize_smiles("CCO")
    assert nm is not None
    events: list[str] = []

    def handler(request: httpx.Request) -> httpx.Response:
        events.append(f"http:{request.method}")
        url = str(request.url)
        if request.method == "GET" and "/inchikey/" in url:
            return httpx.Response(404, text="Not Found")
        if request.method == "POST" and "/smiles/cids/JSON" in url:
            return httpx.Response(200, json={"IdentifierList": {"CID": [702]}})
        if request.method == "GET" and "/categories/compound/702/JSON" in url:
            return httpx.Response(
                200,
                json={
                    "SourceCategories": {
                        "Categories": [
                            {
                                "Category": "Chemical Vendors",
                                "Sources": [{"SourceName": "EthanolVendor", "RegistryID": "ET-702"}],
                            }
                        ]
                    }
                },
            )
        return httpx.Response(404)

    with ResponseCache(tmp_path / "post_cache.sqlite") as cache:
        backend = PubChemVendorBackend(
            cache,
            min_interval_sec=0.0,
            backoff_base_sec=0.0,
            client=httpx.Client(transport=httpx.MockTransport(handler)),
        )
        orig_throttle = backend._throttle

        def tracked_throttle() -> None:
            events.append("throttle")
            orig_throttle()

        monkeypatch.setattr(backend, "_throttle", tracked_throttle)
        offers = backend.lookup(nm.inchikey, nm.canonical_smiles)
        assert len(offers) == 1
        assert offers[0].vendor == "EthanolVendor"
        # Every single outbound HTTP call (GET, POST, GET) must be immediately preceded by _throttle()
        assert events == [
            "throttle",
            "http:GET",
            "throttle",
            "http:POST",
            "throttle",
            "http:GET",
        ]


def test_response_cache_ttl_and_refresh(tmp_path: Path) -> None:
    db = tmp_path / "ttl.sqlite"
    with ResponseCache(db) as cache:
        cache.set("k_old", "pubchem_pug", {"v": 1}, fetched_utc="2020-01-01T00:00:00Z")
        # Infinite TTL (offline mode default) returns old entry
        assert cache.get("k_old", max_age_seconds=None) is not None
        # 1-day TTL treats 2020 entry as expired miss
        assert cache.get("k_old", max_age_seconds=86400.0) is None

    # refresh=True bypasses existing entry until re-populated in current session
    with ResponseCache(db, refresh=True) as refreshing_cache:
        assert refreshing_cache.get("k_old") is None
        refreshing_cache.set("k_old", "pubchem_pug", {"v": 2})
        got = refreshing_cache.get("k_old")
        assert got is not None and got[0] == {"v": 2}


def _mp_cache_writer(args: tuple[str, int, int]) -> int:
    db_path_str, worker_id, n_writes = args
    with ResponseCache(Path(db_path_str)) as c:
        for i in range(n_writes):
            key = f"mp:w{worker_id}:k{i}"
            c.set(key, "pubchem_pug", {"worker": worker_id, "idx": i})
            got = c.get(key)
            assert got is not None and got[0]["idx"] == i
    return n_writes


def test_sqlite_wal_busy_timeout_multithread_multiprocess_and_readonly_uris(
    tmp_path: Path,
) -> None:
    import json
    import multiprocessing as mp
    import os
    import sqlite3
    import stat
    from concurrent.futures import ThreadPoolExecutor

    from buyormake.routes import ReactionLookupBackend
    from buyormake.stock import LocalStockBackend

    cache_db = tmp_path / "shared_cache.sqlite"
    with ResponseCache(cache_db) as cache:
        assert cache._conn is not None
        jmode = cache._conn.execute("PRAGMA journal_mode").fetchone()[0]
        btimeout = cache._conn.execute("PRAGMA busy_timeout").fetchone()[0]
        assert str(jmode).lower() == "wal"
        assert int(btimeout) == 30000

        # 1. >= 8 concurrent threads reading & writing shared ResponseCache
        def _thread_worker(tid: int) -> int:
            for idx in range(25):
                k = f"thr:{tid}:{idx}"
                cache.set(k, "pubchem_pug", {"tid": tid, "idx": idx})
                res = cache.get(k)
                assert res is not None and res[0]["tid"] == tid
            return 25

        with ThreadPoolExecutor(max_workers=10) as pool:
            counts = list(pool.map(_thread_worker, range(10)))
        assert sum(counts) == 250

    # 2. Multi-process concurrent writers against shared cache.sqlite
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=4) as proc_pool:
        mp_counts = proc_pool.map(
            _mp_cache_writer,
            [(str(cache_db), wid, 15) for wid in range(4)],
        )
    assert sum(mp_counts) == 60

    # 3. Read-only mode on ResponseCache skips DDL & writes cleanly
    missing_ro_db = tmp_path / "nonexistent_subdir" / "ro_cache.sqlite"
    with ResponseCache(missing_ro_db, readonly=True) as ro_missing:
        assert ro_missing.get("any") is None
        assert not missing_ro_db.exists()

    with ResponseCache(cache_db, readonly=True) as ro_existing:
        got_mp = ro_existing.get("mp:w0:k0")
        assert got_mp is not None and got_mp[0]["worker"] == 0

    # 4. Read-only stock.sqlite and reactions.sqlite queried across 10 worker threads
    stock_db = tmp_path / "ro_stock.sqlite"
    sc = sqlite3.connect(str(stock_db))
    sc.execute(
        "CREATE TABLE stock (inchikey TEXT, vendor TEXT, source_id TEXT, source_url TEXT, fetched_utc TEXT, license TEXT)"
    )
    eth_ik = normalize_smiles("CCO").inchikey  # type: ignore[union-attr]
    asp_ik = normalize_smiles("CC(=O)Oc1ccccc1C(=O)O").inchikey  # type: ignore[union-attr]
    sc.execute(
        "INSERT INTO stock VALUES (?, 'V1', 'local_stock', 'https://ex/1', '2026-11-18T00:00:00Z', 'CC0-1.0')",
        (eth_ik,),
    )
    sc.commit()
    sc.close()

    rxn_db = tmp_path / "ro_reactions.sqlite"
    rc = sqlite3.connect(str(rxn_db))
    rc.execute(
        "CREATE TABLE reactions (id INTEGER PRIMARY KEY, product_inchikey TEXT, product_smiles TEXT, "
        "reactants_json TEXT, reaction_smiles TEXT, reference TEXT, source TEXT)"
    )
    rc.execute(
        "INSERT INTO reactions (product_inchikey, product_smiles, reactants_json, reaction_smiles, reference, source) "
        "VALUES (?, 'CC(=O)Oc1ccccc1C(=O)O', ?, 'CC(=O)O.O=C(O)c1ccccc1O>>CC(=O)Oc1ccccc1C(=O)O', 'ref1', 'uspto_lookup')",
        (asp_ik, json.dumps(["CC(=O)O", "O=C(O)c1ccccc1O"])),
    )
    rc.commit()
    rc.close()

    os.chmod(stock_db, stat.S_IREAD)
    os.chmod(rxn_db, stat.S_IREAD)
    try:
        stock_backend = LocalStockBackend([stock_db])
        rxn_backend = ReactionLookupBackend(rxn_db)
        try:
            with ThreadPoolExecutor(max_workers=10) as pool:
                s_res = list(pool.map(lambda _: len(stock_backend.lookup(eth_ik)), range(30)))
                r_res = list(
                    pool.map(
                        lambda _: len(rxn_backend.one_step(asp_ik, "CC(=O)Oc1ccccc1C(=O)O")),
                        range(30),
                    )
                )
            assert all(n == 1 for n in s_res)
            assert all(n == 1 for n in r_res)
        finally:
            stock_backend.close()
            rxn_backend.close()
    finally:
        os.chmod(stock_db, stat.S_IWRITE | stat.S_IREAD)
        os.chmod(rxn_db, stat.S_IWRITE | stat.S_IREAD)
