"""Tests for PubChem caching and offline mode (mocked HTTP; zero network calls)."""

from __future__ import annotations

import json
from pathlib import Path

import httpx
import pytest

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
