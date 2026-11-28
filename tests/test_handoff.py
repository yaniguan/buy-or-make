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
                est_cost=Quantity(value=1.25, unit="USD/g"),
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
                est_cost=Quantity(value=4.2, unit="USD/g"),
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


def test_evaluate_candidates_persists_leaf_offers_for_automatic_handoff(tmp_path) -> None:
    from buyormake import build_precursor_attribution_basis, evaluate_candidates

    stock_tsv = tmp_path / "stock.tsv"
    stock_tsv.write_text(
        "smiles\tvendor\tcatalog_id\tpack_value\tpack_unit\tprice\tcurrency\tlead_time_days\tin_stock\tsource_id\tsource_url\tlicense\n"
        "COc1ccc(C(=O)O)cc1\tRefV\tMBA-1\t1.0\tg\t0.40\tUSD\t2\ttrue\tlocal_stock\thttps://example.org/mba\tCC0-1.0\n"
        "Nc1ccc(Br)cc1\tRefV\tBRA-2\t1.0\tg\t0.80\tUSD\t2\ttrue\tlocal_stock\thttps://example.org/bra\tCC0-1.0\n"
        "OB(O)c1ccccn1\tRefV\tPYB-3\t1.0\tg\t3.50\tUSD\t3\ttrue\tlocal_stock\thttps://example.org/pyb\tCC0-1.0\n",
        encoding="utf-8",
    )

    _, records = evaluate_candidates(
        ["COc1ccc(C(=O)Nc2ccc(-c3ccccn3)cc2)cc1"],
        candidate_names=["Analog-3"],
        max_steps=2,
        backends=("local_stock",),
        route_backend="template",
        offline=True,
        cache_path=tmp_path / "cache.sqlite",
        stock_paths=[stock_tsv],
        reaction_db_path=tmp_path / "no_rxn.sqlite",
        scscore_weights=tmp_path / "no_sc.json.gz",
    )
    rec = records[0]
    assert len(rec.leaf_offers) == 3
    assert rec.lead_time_days == 3
    assert rec.moq_pack_g == pytest.approx(1.0)
    assert len(rec.precursor_smiles) == 3

    # Calling to_procurement_constraint(rec) WITHOUT leaf_offers_by_inchikey uses rec.leaf_offers!
    c_auto = to_procurement_constraint(rec)
    assert c_auto.schedulable is True
    assert c_auto.lead_time_days == 3
    assert c_auto.lead_time_complete is True
    assert c_auto.release_time_minutes == 4800
    assert c_auto.moq_pack_g == pytest.approx(1.0)
    assert c_auto.min_order_outlay_usd == pytest.approx(4.70)
    assert build_labshift_release_map([rec]) == {"Analog-3": 4800}

    basis, rxn_steps = build_precursor_attribution_basis(rec.routes[0], rec.leaf_offers)
    assert {b["catalog_id"] for b in basis} == {"MBA-1", "BRA-2", "PYB-3"}
    assert len(rxn_steps) == 2


def test_handoff_respects_in_stock_false_and_preserves_in_stock_none() -> None:
    """Issue #36: all-out-of-stock offers give schedulable=False, release_time_minutes=None, omitted from release map."""
    asp_nm = normalize_smiles("CC(=O)Oc1ccccc1C(=O)O")
    etoh_nm = normalize_smiles("CCO")
    assert asp_nm and etoh_nm

    rec_oos = AvailabilityRecord(
        molecule=Molecule(smiles=asp_nm.canonical_smiles, inchikey=asp_nm.inchikey),
        offers=[
            Offer(
                inchikey=asp_nm.inchikey,
                vendor="OutVendor",
                pack_size=Quantity(value=1.0, unit="g", raw_unit_str="1 g"),
                price=999.0,
                currency="USD",
                price_basis="listed_catalog",
                in_stock=False,
                source_id="local_stock",
                source_url="https://example.org/oos",
                license="CC0-1.0",
            )
        ],
        licenses=["CC0-1.0"],
    ).sealed()

    c_oos = to_procurement_constraint(rec_oos)
    assert c_oos.status == "unknown"
    assert c_oos.schedulable is False
    assert c_oos.release_time_minutes is None
    assert c_oos.est_cost is not None and c_oos.est_cost.value == pytest.approx(999.0)
    assert any("flagged out of stock" in note for note in c_oos.notes)
    assert not any("unlisted" in note for note in c_oos.notes)
    assert build_labshift_release_map([rec_oos]) == {}

    rec_none_flag = AvailabilityRecord(
        molecule=Molecule(smiles=etoh_nm.canonical_smiles, inchikey=etoh_nm.inchikey, name="Ethanol"),
        offers=[
            Offer(
                inchikey=etoh_nm.inchikey,
                vendor="PubVendor",
                lead_time_days=3,
                in_stock=None,
                source_id="pubchem_pug",
                source_url="https://example.org/pub",
                license="custom/terms",
            )
        ],
        licenses=["custom/terms"],
    ).sealed()
    c_none = to_procurement_constraint(rec_none_flag)
    assert c_none.status == "buyable"
    assert c_none.schedulable is True
    assert c_none.release_time_minutes == 3 * 1440
    assert build_labshift_release_map([rec_oos, rec_none_flag]) == {"Ethanol": 3 * 1440}
