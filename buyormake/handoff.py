"""Handoff adapters connecting `buy-or-make` (`L5`) to `assayfold` (`L0`), `labshift` (`L5`), and `lab-decide`."""

from __future__ import annotations

from typing import Iterable, Literal, Optional

from pydantic import Field

from .canonical import normalize_smiles
from .schema import (
    AssayfoldRecord,
    AvailabilityRecord,
    Offer,
    Quantity,
    _Model,
    select_record_license,
)


class ProcurementScheduleConstraint(_Model):
    """Scheduling release-time and cost handoff record for `labshift` / `lab-decide`."""

    molecule_id: str
    smiles: str
    inchikey: Optional[str] = None
    status: Literal["buyable", "makeable", "unknown"]
    schedulable: bool
    lead_time_days: Optional[int] = None
    n_synthesis_steps: int = 0
    release_time_minutes: Optional[int] = None
    lead_time_complete: bool = False
    est_cost: Optional[Quantity] = None
    governing_license: str = "MIT"
    notes: list[str] = Field(default_factory=list)


def _min_offer_lead_days(offers: Iterable[Offer]) -> Optional[int]:
    """Return shortest reported lead_time_days among offers not explicitly out of stock."""
    in_stock_leads = [
        off.lead_time_days
        for off in offers
        if off.lead_time_days is not None and off.in_stock is not False
    ]
    if in_stock_leads:
        return min(in_stock_leads)
    any_leads = [off.lead_time_days for off in offers if off.lead_time_days is not None]
    return min(any_leads) if any_leads else None


def to_procurement_constraint(
    record: AvailabilityRecord,
    *,
    leaf_offers_by_inchikey: Optional[dict[str, list[Offer]]] = None,
    minutes_per_day: int = 1440,
    minutes_per_synthesis_step: int = 240,
    default_lead_time_days: int = 0,
) -> ProcurementScheduleConstraint:
    """Map an `AvailabilityRecord` (`Offer.lead_time_days` & `Route.n_steps`) to a `labshift` `release_time` offset.

    Semantics:
    - `buyable`: earliest experiment `release_time` is `min(Offer.lead_time_days) * minutes_per_day`
      (falling back to `default_lead_time_days` with `lead_time_complete=False` if no offer lists lead time).
    - `makeable`: synthesis cannot start until *all* leaf building blocks arrive (`max` over leaves of
      each leaf's minimum vendor `lead_time_days`), plus `n_steps * minutes_per_synthesis_step`.
    - `unknown`: `schedulable=False` and `release_time_minutes=None`.
    """
    if minutes_per_day <= 0:
        raise ValueError("minutes_per_day must be > 0")
    if minutes_per_synthesis_step < 0:
        raise ValueError("minutes_per_synthesis_step must be >= 0")
    if default_lead_time_days < 0:
        raise ValueError("default_lead_time_days must be >= 0")

    from .cost import cheapest_offer_for_target

    mol_id = record.molecule.name or record.molecule.inchikey or record.molecule.smiles
    gov_lic = select_record_license(
        [*record.licenses, *(o.license for o in record.offers if o.license)],
        default="MIT",
    )
    notes: list[str] = []

    if record.offers:
        lead = _min_offer_lead_days(record.offers)
        complete = lead is not None
        eff_lead = lead if lead is not None else default_lead_time_days
        if not complete:
            notes.append(
                f"vendor lead_time_days unlisted; using default_lead_time_days={default_lead_time_days}"
            )
        release_min = eff_lead * minutes_per_day
        return ProcurementScheduleConstraint(
            molecule_id=mol_id,
            smiles=record.molecule.smiles,
            inchikey=record.molecule.inchikey,
            status="buyable",
            schedulable=True,
            lead_time_days=lead,
            n_synthesis_steps=0,
            release_time_minutes=release_min,
            lead_time_complete=complete,
            est_cost=cheapest_offer_for_target(record.offers),
            governing_license=gov_lic,
            notes=notes,
        )

    solved = [r for r in record.routes if r.leaves and all(r.leaves_purchasable)]
    if solved:
        best_route = solved[0]
        leaf_map = leaf_offers_by_inchikey or {}
        leaf_leads: list[int] = []
        missing_lead_leaves: list[str] = []
        for leaf_smi in dict.fromkeys(best_route.leaves):
            nm = normalize_smiles(leaf_smi)
            offs = leaf_map.get(nm.inchikey, []) if nm is not None else []
            l_days = _min_offer_lead_days(offs)
            if l_days is None:
                missing_lead_leaves.append(leaf_smi)
            else:
                leaf_leads.append(l_days)

        complete = len(missing_lead_leaves) == 0 and len(best_route.leaves) > 0
        bottleneck_days = max(leaf_leads) if leaf_leads else default_lead_time_days
        if missing_lead_leaves:
            if leaf_leads:
                bottleneck_days = max(bottleneck_days, default_lead_time_days)
            notes.append(
                f"{len(missing_lead_leaves)} leaf/leaves lack lead_time_days; "
                f"applied default_lead_time_days={default_lead_time_days} floor"
            )
        synth_min = best_route.n_steps * minutes_per_synthesis_step
        release_min = bottleneck_days * minutes_per_day + synth_min
        return ProcurementScheduleConstraint(
            molecule_id=mol_id,
            smiles=record.molecule.smiles,
            inchikey=record.molecule.inchikey,
            status="makeable",
            schedulable=True,
            lead_time_days=max(leaf_leads) if complete else None,
            n_synthesis_steps=best_route.n_steps,
            release_time_minutes=release_min,
            lead_time_complete=complete,
            est_cost=best_route.est_cost,
            governing_license=gov_lic,
            notes=notes,
        )

    notes.append("unschedulable: molecule is neither in stock nor makeable within search bounds")
    return ProcurementScheduleConstraint(
        molecule_id=mol_id,
        smiles=record.molecule.smiles,
        inchikey=record.molecule.inchikey,
        status="unknown",
        schedulable=False,
        lead_time_days=None,
        n_synthesis_steps=0,
        release_time_minutes=None,
        lead_time_complete=False,
        est_cost=None,
        governing_license=gov_lic,
        notes=notes,
    )


def build_labshift_release_map(
    records: Iterable[AvailabilityRecord],
    *,
    leaf_offers_by_inchikey: Optional[dict[str, list[Offer]]] = None,
    minutes_per_day: int = 1440,
    minutes_per_synthesis_step: int = 240,
    default_lead_time_days: int = 0,
) -> dict[str, int]:
    """Build `{molecule_id: release_time_minutes}` map for schedulable candidates in `labshift`."""
    out: dict[str, int] = {}
    for rec in records:
        c = to_procurement_constraint(
            rec,
            leaf_offers_by_inchikey=leaf_offers_by_inchikey,
            minutes_per_day=minutes_per_day,
            minutes_per_synthesis_step=minutes_per_synthesis_step,
            default_lead_time_days=default_lead_time_days,
        )
        if c.schedulable and c.release_time_minutes is not None:
            out[c.molecule_id] = c.release_time_minutes
    return out


__all__ = [
    "AssayfoldRecord",
    "ProcurementScheduleConstraint",
    "build_labshift_release_map",
    "select_record_license",
    "to_procurement_constraint",
]
