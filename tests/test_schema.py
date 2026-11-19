"""Tests for schema validation, extra='forbid', deterministic hashing, and source registry."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pydantic import ValidationError

from buyormake import (
    SOURCES,
    AvailabilityRecord,
    Molecule,
    Offer,
    Provenance,
    Quantity,
    RankResult,
    Route,
    RouteStep,
    export_json_schema,
)


def test_extra_forbid_rejects_unknown_fields() -> None:
    with pytest.raises(ValidationError):
        Molecule(smiles="CCO", extra_field="nope")  # type: ignore[call-arg]

    with pytest.raises(ValidationError):
        Offer(
            inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            vendor="TestVendor",
            source_id="local_stock",
            source_url="https://example.org",
            license="CC0-1.0",
            invented_margin=0.5,  # type: ignore[call-arg]
        )

    with pytest.raises(ValidationError):
        RankResult(candidate="CCO", status="buyable", bogus=1)  # type: ignore[call-arg]


def test_record_hash_deterministic_and_excludes_recorded_utc() -> None:
    mol = Molecule(smiles="CCO", inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N", formula="C2H6O")
    rec1 = AvailabilityRecord(
        molecule=mol,
        provenance=Provenance(recorded_utc="2026-01-01T00:00:00Z"),
    ).sealed()
    rec2 = AvailabilityRecord(
        molecule=mol,
        provenance=Provenance(recorded_utc="2026-11-18T12:34:56Z"),
    ).sealed()

    assert rec1.record_hash is not None
    assert len(rec1.record_hash) == 64
    assert rec1.record_hash == rec2.record_hash

    rec3 = AvailabilityRecord(
        molecule=mol,
        offers=[
            Offer(
                inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                vendor="VendorA",
                price=None,
                currency=None,
                price_basis="unknown",
                source_id="local_stock",
                source_url="https://example.org",
                fetched_utc="2026-11-18T00:00:00Z",
                license="CC0-1.0",
            )
        ],
    ).sealed()
    assert rec3.record_hash != rec1.record_hash


def test_all_sources_valid() -> None:
    assert len(SOURCES) >= 10
    for src in SOURCES:
        assert src.license_evidence.startswith("http")
        assert src.license_checked
        assert src.status in {"integrated", "fetch_only", "not_used", "planned"}


def test_zinc_indirect_license_check_is_machine_visible() -> None:
    by_id = {src.id: src for src in SOURCES}
    zinc = by_id["zinc"]
    assert "indirect" in zinc.license_checked
    assert zinc.notes is not None and "LICENSE CHECK INDIRECT" in zinc.notes
    assert all("indirect" not in src.license_checked for src in SOURCES if src.id != "zinc")


def test_export_json_schema(tmp_path: Path) -> None:
    paths = export_json_schema(tmp_path)
    assert len(paths) == 5
    for p in paths:
        data = json.loads(p.read_text(encoding="utf-8"))
        assert "properties" in data
        assert data.get("additionalProperties") is False


def test_route_and_quantity_models() -> None:
    r = Route(
        target_inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        n_steps=1,
        steps=[
            RouteStep(
                product="CC(=O)Oc1ccccc1C(=O)O",
                reactants=["CC(=O)O", "O=C(O)c1ccccc1O"],
                reaction_smiles="CC(=O)O.O=C(O)c1ccccc1O>>CC(=O)Oc1ccccc1C(=O)O",
                template_or_reference="esterification",
                source="template_fallback",
            )
        ],
        leaves=["CC(=O)O", "O=C(O)c1ccccc1O"],
        leaves_purchasable=[True, True],
        method="template_fallback",
        est_cost=Quantity(value=0.12, unit="USD"),
        est_cost_method="leaf sum",
        est_cost_complete=True,
    )
    assert r.n_steps == 1
