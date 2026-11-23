"""Tests for schema validation, extra='forbid', cross-field invariants, zero drift, and source registry."""

from __future__ import annotations

import json
import re
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
    Source,
    export_json_schema,
)
from buyormake.schema import Uncertainty


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
    iso_date_re = re.compile(r"^\d{4}-\d{2}-\d{2}$")
    for src in SOURCES:
        assert src.license_evidence.startswith("http")
        assert iso_date_re.match(src.license_checked)
        assert src.license_check_method in {"direct", "indirect"}
        assert src.status in {"integrated", "fetch_only", "not_used", "planned"}


def test_zinc_indirect_license_check_is_machine_visible() -> None:
    by_id = {src.id: src for src in SOURCES}
    zinc = by_id["zinc"]
    assert zinc.license_checked == "2026-11-18"
    assert zinc.license_check_method == "indirect"
    assert zinc.notes is not None and "LICENSE CHECK INDIRECT" in zinc.notes
    assert all(src.license_check_method == "direct" for src in SOURCES if src.id != "zinc")

    with pytest.raises(ValidationError):
        Source(
            id="bad_date",
            name="Bad Date",
            homepage="https://example.org",
            license="MIT",
            license_evidence="https://example.org/license",
            license_checked="2026-11-18 (free text not allowed)",
            commercial_use_ok=True,
            redistribution_ok=True,
            attribution="Test",
            access="github",
        )


def test_offer_and_quantity_cross_field_invariants() -> None:
    base = {
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        "vendor": "TestVendor",
        "source_id": "local_stock",
        "source_url": "https://example.org",
        "license": "CC0-1.0",
    }
    # Negative price rejected
    with pytest.raises(ValidationError):
        Offer(**base, price=-5.0, currency="USD", price_basis="listed_catalog")
    # Price without currency rejected
    with pytest.raises(ValidationError):
        Offer(**base, price=10.0, currency=None, price_basis="listed_catalog")
    # Currency without price rejected
    with pytest.raises(ValidationError):
        Offer(**base, price=None, currency="USD", price_basis="unknown")
    # listed_catalog with price=None rejected
    with pytest.raises(ValidationError):
        Offer(**base, price=None, currency=None, price_basis="listed_catalog")
    # Negative lead_time_days rejected
    with pytest.raises(ValidationError):
        Offer(**base, lead_time_days=-3)
    # Non-positive pack_size rejected
    with pytest.raises(ValidationError):
        Offer(**base, pack_size=Quantity(value=0.0, unit="g"))
    # NaN / Inf in Quantity or Offer rejected
    with pytest.raises(ValidationError):
        Quantity(value=float("nan"), unit="USD")
    with pytest.raises(ValidationError):
        Quantity(value=float("inf"), unit="USD")
    with pytest.raises(ValidationError):
        Offer(**base, price=float("nan"), currency="USD", price_basis="listed_catalog")
    # Uncertainty low > high rejected
    with pytest.raises(ValidationError):
        Uncertainty(kind="interval", low=5.0, high=2.0)


def test_export_json_schema_matches_committed_files(tmp_path: Path) -> None:
    paths = export_json_schema(tmp_path)
    assert len(paths) == 5
    committed_dir = Path(__file__).resolve().parents[1] / "schema"
    for p in paths:
        generated = p.read_text(encoding="utf-8")
        data = json.loads(generated)
        assert "properties" in data
        assert data.get("additionalProperties") is False
        committed_file = committed_dir / p.name
        assert committed_file.exists(), f"Missing committed schema file {committed_file}"
        assert committed_file.read_text(encoding="utf-8") == generated, (
            f"Schema drift detected in {p.name}; run `buy-or-make export-schema --out schema`"
        )


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
    with pytest.raises(ValidationError):
        Route(
            target_inchikey="BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            n_steps=1,
            steps=[],
            leaves=["CC(=O)O"],
            leaves_purchasable=[True, False],
            method="template_fallback",
            est_cost_method="leaf sum",
            est_cost_complete=False,
        )
