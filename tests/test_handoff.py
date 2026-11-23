"""Tests for assayfold Record conversion and labshift / lab-decide procurement handoff."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from buyormake import (
    AssayfoldRecord,
    AvailabilityRecord,
    Molecule,
    Offer,
    Provenance,
    Quantity,
    Route,
    RouteStep,
    SAScores,
    build_labshift_release_map,
    select_record_license,
    to_procurement_constraint,
)
from buyormake.canonical import normalize_smiles


def test_to_assayfold_record_roundtrip_and_license_gate() -> None:
    # Restrictiveness ordering: custom/terms > CC-BY-SA-4.0 > CC-BY-4.0 > BSD-3-Clause > MIT > CC0-1.0
    assert select_record_license(["MIT", "BSD-3-Clause", "CC-BY-SA-4.0"]) == "CC-BY-SA-4.0"
    assert select_record_license(["CC-BY-SA-4.0", "custom/terms"]) == "custom/terms"

    nm = normalize_smiles("CC(=O)Nc1ccc(F)cc1")
    assert nm is not None
    avail = AvailabilityRecord(
        molecule=Molecule(
            smiles=nm.canonical_smiles,
            inchikey=nm.inchikey,
            name="4-Fluoroacetanilide",
            formula=nm.formula,
        ),
        routes=[
            Route(
                target_inchikey=nm.inchikey,
                n_steps=1,
                steps=[
                    RouteStep(
                        product=nm.canonical_smiles,
                        reactants=["CC(=O)O", "Nc1ccc(F)cc1"],
                        template_or_reference="ord-001",
                        source="ord_lookup",
                    )
                ],
                leaves=["CC(=O)O", "Nc1ccc(F)cc1"],
                leaves_purchasable=[True, True],
                method="reaction_lookup (ord_lookup)",
                score=0.88,
                est_cost=Quantity(value=1.25, unit="USD"),
                est_cost_method="leaf sum",
                est_cost_complete=True,
            )
        ],
        sa_scores=SAScores(sascore=1.42, scscore=1.65),
        provenance=Provenance(recorded_utc="2026-11-18T00:00:00Z", used=["ord_data"]),
        licenses=["BSD-3-Clause", "CC-BY-SA-4.0", "MIT"],
    ).sealed()

    af_rec = avail.to_assayfold_record()
    assert af_rec.record_id.startswith("af:")
    assert len(af_rec.record_id) == 27
    assert af_rec.kind == "molecule"
    assert af_rec.payload.smiles == nm.canonical_smiles
    assert af_rec.license == "CC-BY-SA-4.0"
    assert af_rec.properties["min_synthesis_steps"].value == 1.0
    assert af_rec.properties["est_route_cost"].value == pytest.approx(1.25)
    assert af_rec.properties["route_score"].value == pytest.approx(0.88)
    assert af_rec.properties["sascore"].value == pytest.approx(1.42)
    assert avail.record_hash in af_rec.provenance.derived_from
    assert "status:makeable" in af_rec.tags

    # Round-trip JSON validation & timestamp-invariant content hash
    reloaded = AssayfoldRecord.model_validate_json(af_rec.model_dump_json())
    assert reloaded.record_hash == af_rec.record_hash == reloaded.content_hash()
    reloaded_later = reloaded.model_copy(
        update={"provenance": reloaded.provenance.model_copy(update={"recorded_utc": "2026-12-31T23:59:59Z"})}
    )
    assert reloaded_later.content_hash() == af_rec.record_hash

    with pytest.raises(ValidationError):
        AssayfoldRecord.model_validate({**af_rec.model_dump(mode="json"), "extra_key": 1})


def test_labshift_procurement_release_time_handoff() -> None:
    # 1. Buyable candidate with 2-day vs 5-day vendor offers -> earliest release = 2 * 1440 = 2880 min
    aspirin_nm = normalize_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert aspirin_nm is not None
    rec_buy = AvailabilityRecord(
        molecule=Molecule(smiles=aspirin_nm.canonical_smiles, inchikey=aspirin_nm.inchikey, name="Aspirin"),
        offers=[
            Offer(
                inchikey=aspirin_nm.inchikey,
                vendor="SlowVendor",
                lead_time_days=5,
                in_stock=True,
                source_id="local_stock",
                source_url="https://example.org/slow",
                license="CC0-1.0",
            ),
            Offer(
                inchikey=aspirin_nm.inchikey,
                vendor="FastVendor",
                lead_time_days=2,
                in_stock=True,
                source_id="local_stock",
                source_url="https://example.org/fast",
                license="CC0-1.0",
            ),
        ],
        licenses=["CC0-1.0"],
    ).sealed()

    # 2. Makeable 2-step candidate with leaf lead times 1 day and 4 days -> bottleneck = 4 days + 2 * 240 min
    target_nm = normalize_smiles("COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1")
    l1_nm = normalize_smiles("COc1ccc(C(=O)O)cc1")
    l2_nm = normalize_smiles("Nc1ccc(Br)cc1")
    l3_nm = normalize_smiles("OB(O)c1ccccn1")
    assert target_nm and l1_nm and l2_nm and l3_nm

    rec_make = AvailabilityRecord(
        molecule=Molecule(smiles=target_nm.canonical_smiles, inchikey=target_nm.inchikey, name="Analog-3"),
        routes=[
            Route(
                target_inchikey=target_nm.inchikey,
                n_steps=2,
                steps=[
                    RouteStep(
                        product=target_nm.canonical_smiles,
                        reactants=[l1_nm.canonical_smiles, "Nc1ccc(-c2ccccn2)cc1"],
                        template_or_reference="amide",
                        source="template_fallback",
                    ),
                    RouteStep(
                        product="Nc1ccc(-c2ccccn2)cc1",
                        reactants=[l2_nm.canonical_smiles, l3_nm.canonical_smiles],
                        template_or_reference="suzuki",
                        source="template_fallback",
                    ),
                ],
                leaves=[l1_nm.canonical_smiles, l2_nm.canonical_smiles, l3_nm.canonical_smiles],
                leaves_purchasable=[True, True, True],
                method="template_fallback",
                est_cost=Quantity(value=4.2, unit="USD"),
                est_cost_method="leaf sum",
                est_cost_complete=True,
            )
        ],
        licenses=["CC0-1.0"],
    ).sealed()

    # 3. Unknown candidate -> unschedulable (excluded from release map)
    cubane_nm = normalize_smiles("C12C3C4C1C5C2C3C45")
    assert cubane_nm is not None
    rec_unk = AvailabilityRecord(
        molecule=Molecule(smiles=cubane_nm.canonical_smiles, inchikey=cubane_nm.inchikey, name="Cubane"),
    ).sealed()

    leaf_offers = {
        l1_nm.inchikey: [
            Offer(
                inchikey=l1_nm.inchikey,
                vendor="V1",
                lead_time_days=1,
                in_stock=True,
                source_id="local_stock",
                source_url="https://example.org/l1",
                license="CC0-1.0",
            )
        ],
        l2_nm.inchikey: [
            Offer(
                inchikey=l2_nm.inchikey,
                vendor="V2",
                lead_time_days=4,
                in_stock=True,
                source_id="local_stock",
                source_url="https://example.org/l2",
                license="CC0-1.0",
            )
        ],
        l3_nm.inchikey: [
            Offer(
                inchikey=l3_nm.inchikey,
                vendor="V3",
                lead_time_days=2,
                in_stock=True,
                source_id="local_stock",
                source_url="https://example.org/l3",
                license="CC0-1.0",
            )
        ],
    }

    c_make = to_procurement_constraint(
        rec_make,
        leaf_offers_by_inchikey=leaf_offers,
        minutes_per_day=1440,
        minutes_per_synthesis_step=240,
    )
    assert c_make.schedulable is True
    assert c_make.lead_time_complete is True
    assert c_make.lead_time_days == 4
    assert c_make.release_time_minutes == 4 * 1440 + 2 * 240

    c_unk = to_procurement_constraint(rec_unk)
    assert c_unk.schedulable is False
    assert c_unk.release_time_minutes is None

    rel_map = build_labshift_release_map(
        [rec_buy, rec_make, rec_unk],
        leaf_offers_by_inchikey=leaf_offers,
        minutes_per_day=1440,
        minutes_per_synthesis_step=240,
    )
    assert rel_map == {
        "Aspirin": 2 * 1440,
        "Analog-3": 4 * 1440 + 2 * 240,
    }
