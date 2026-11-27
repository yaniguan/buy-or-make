"""Handoff adapters connecting `buy-or-make` (`L5`) to `assayfold` (`L0`), `labshift` (`L5`), `lab-attrib` (`L4`), and `lab-decide`."""

from __future__ import annotations

from collections import Counter
from typing import Any, Iterable, Literal, Mapping, Optional

from pydantic import Field
from rdkit import Chem

from .canonical import normalize_smiles
from .cost import (
    cheapest_offer_for_target,
    estimate_route_cost_bounds,
    min_outlay_for_offers,
)
from .retro_star import RETRO_STAR_TEMPLATE_DELTAS
from .routes import RETRO_TEMPLATE_DELTAS
from .schema import (
    AssayfoldRecord,
    AvailabilityRecord,
    Offer,
    Quantity,
    Route,
    _Model,
    select_record_license,
)

_ROLE_BORONIC = Chem.MolFromSmarts("[#6]B([OH])[OH]")
_ROLE_SULFONYL_CL = Chem.MolFromSmarts("[#6]S(=O)(=O)Cl")
_ROLE_ISOCYANATE = Chem.MolFromSmarts("[#6]N=C=O")
_ROLE_CARBOXYLIC = Chem.MolFromSmarts("[CX3](=O)[OH1]")
_ROLE_ARYL_HALIDE = Chem.MolFromSmarts("[c,CX3][F,Cl,Br,I]")
_ROLE_ALKYL_HALIDE = Chem.MolFromSmarts("[CX4][Cl,Br,I]")
_ROLE_ALDEHYDE = Chem.MolFromSmarts("[CX3H1](=O)[#6,#1]")
_ROLE_AMINE = Chem.MolFromSmarts("[N;!a;H2,H1;!$(NC(=[O,S]));!$(NS(=O)=O)]")
_ROLE_ALCOHOL = Chem.MolFromSmarts("[OH1;!$(OC(=[O,S]));!$(OB);!$(OS)]")


def _infer_synthon_role(mol: Optional[Chem.Mol]) -> str:
    if mol is None:
        return "building_block"
    if mol.HasSubstructMatch(_ROLE_BORONIC):
        return "boronic_acid"
    if mol.HasSubstructMatch(_ROLE_SULFONYL_CL):
        return "sulfonyl_chloride"
    if mol.HasSubstructMatch(_ROLE_ISOCYANATE):
        return "isocyanate"
    if mol.HasSubstructMatch(_ROLE_CARBOXYLIC):
        return "carboxylic_acid"
    if mol.HasSubstructMatch(_ROLE_ALDEHYDE):
        return "aldehyde"
    if mol.HasSubstructMatch(_ROLE_ARYL_HALIDE):
        return "aryl_halide"
    if mol.HasSubstructMatch(_ROLE_ALKYL_HALIDE):
        return "alkyl_halide"
    if mol.HasSubstructMatch(_ROLE_AMINE):
        return "amine_nucleophile"
    if mol.HasSubstructMatch(_ROLE_ALCOHOL):
        return "alcohol_nucleophile"
    return "building_block"


def _heavy_counts_from_smiles(smiles: str) -> Counter[str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return Counter()
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


class ProcurementScheduleConstraint(_Model):
    """Scheduling release-time, MOQ pack, and cost handoff record for `labshift` / `lab-decide`."""

    molecule_id: str
    smiles: str
    inchikey: Optional[str] = None
    status: Literal["buyable", "makeable", "unknown"]
    schedulable: bool
    lead_time_days: Optional[int] = None
    n_synthesis_steps: int = 0
    release_time_minutes: Optional[int] = None
    lead_time_complete: bool = False
    moq_pack_g: Optional[float] = None
    selected_pack_size: Optional[Quantity] = None
    min_order_outlay_usd: Optional[float] = None
    precursor_smiles: list[str] = Field(default_factory=list)
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


def build_precursor_attribution_basis(
    route: Route,
    leaf_offers_by_inchikey: Mapping[str, list[Offer]],
    *,
    currency: str = "USD",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Build structured `(precursor_basis, reaction_steps)` records for `lab-attrib` (`L4`)."""
    leaf_counts = Counter(route.leaves)
    precursor_basis: list[dict[str, Any]] = []
    for raw_leaf in dict.fromkeys(route.leaves):
        lnm = normalize_smiles(raw_leaf)
        can_smi = lnm.canonical_smiles if lnm is not None else raw_leaf
        ik = lnm.inchikey if lnm is not None else None
        offs = list(leaf_offers_by_inchikey.get(ik, [])) if ik is not None else []
        _, sel_off, moq_g = min_outlay_for_offers(offs, currency=currency, required_grams=1.0)
        precursor_basis.append(
            {
                "smiles": can_smi,
                "inchikey": ik,
                "catalog_id": sel_off.catalog_id if sel_off is not None else None,
                "vendor": sel_off.vendor if sel_off is not None else None,
                "role": _infer_synthon_role(lnm.mol if lnm is not None else None),
                "stoich_equiv": float(leaf_counts[raw_leaf]),
                "lead_time_days": _min_offer_lead_days(offs),
                "moq_pack_g": moq_g,
            }
        )

    reaction_steps: list[dict[str, Any]] = []
    for idx, step in enumerate(route.steps):
        ref = step.template_or_reference or ""
        tpl_name = ref.split(":", 1)[0].strip()
        if tpl_name in RETRO_STAR_TEMPLATE_DELTAS:
            delta = dict(RETRO_STAR_TEMPLATE_DELTAS[tpl_name])
        elif tpl_name in RETRO_TEMPLATE_DELTAS:
            delta = dict(RETRO_TEMPLATE_DELTAS[tpl_name])
        else:
            p_counts = _heavy_counts_from_smiles(step.product)
            r_counts: Counter[str] = Counter()
            for r_smi in step.reactants:
                r_counts.update(_heavy_counts_from_smiles(r_smi))
            delta = {
                sym: int(r_counts[sym] - p_counts[sym])
                for sym in sorted(set(r_counts) | set(p_counts))
                if r_counts[sym] - p_counts[sym] != 0
            }
        reaction_steps.append(
            {
                "step_index": idx,
                "product": step.product,
                "reactants": list(step.reactants),
                "template_or_reference": step.template_or_reference,
                "source": step.source,
                "leaving_delta": delta,
            }
        )
    return precursor_basis, reaction_steps


def to_procurement_constraint(
    record: AvailabilityRecord,
    *,
    leaf_offers_by_inchikey: Optional[dict[str, list[Offer]]] = None,
    minutes_per_day: int = 1440,
    minutes_per_synthesis_step: int = 240,
    default_lead_time_days: int = 0,
    target_scale_g: float = 1.0,
) -> ProcurementScheduleConstraint:
    """Map an `AvailabilityRecord` (`Offer.lead_time_days` & `Route.n_steps`) to a `labshift` `release_time` offset.

    Automatically consumes `record.leaf_offers` for `makeable` routes when `leaf_offers_by_inchikey`
    is omitted, and merges explicit `leaf_offers_by_inchikey` overrides when supplied.
    """
    if minutes_per_day <= 0:
        raise ValueError("minutes_per_day must be > 0")
    if minutes_per_synthesis_step < 0:
        raise ValueError("minutes_per_synthesis_step must be >= 0")
    if default_lead_time_days < 0:
        raise ValueError("default_lead_time_days must be >= 0")
    if target_scale_g <= 0:
        raise ValueError("target_scale_g must be > 0")

    mol_id = record.molecule.name or record.molecule.inchikey or record.molecule.smiles
    gov_lic = select_record_license(
        [
            *record.licenses,
            *(o.license for o in record.offers if o.license),
            *(lo.license for offs in record.leaf_offers.values() for lo in offs if lo.license),
        ],
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
        outlay_q, sel_off, moq_g = min_outlay_for_offers(
            record.offers, currency="USD", required_grams=target_scale_g
        )
        eff_moq_g = moq_g if moq_g is not None else record.moq_pack_g
        sel_pack = sel_off.pack_size if (sel_off is not None and sel_off.pack_size is not None) else None
        outlay_usd = outlay_q.value if outlay_q is not None else None
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
            moq_pack_g=eff_moq_g,
            selected_pack_size=sel_pack,
            min_order_outlay_usd=outlay_usd,
            precursor_smiles=[],
            est_cost=cheapest_offer_for_target(record.offers),
            governing_license=gov_lic,
            notes=notes,
        )

    solved = [r for r in record.routes if r.leaves and all(r.leaves_purchasable)]
    if solved:
        best_route = solved[0]
        leaf_map: dict[str, list[Offer]] = dict(record.leaf_offers)
        if leaf_offers_by_inchikey:
            leaf_map.update(leaf_offers_by_inchikey)

        leaf_leads: list[int] = []
        missing_lead_leaves: list[str] = []
        precursors: list[str] = []
        for leaf_smi in dict.fromkeys(best_route.leaves):
            nm = normalize_smiles(leaf_smi)
            precursors.append(nm.canonical_smiles if nm is not None else leaf_smi)
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

        breakdown = estimate_route_cost_bounds(
            record.molecule.smiles,
            best_route.leaves,
            leaf_map,
            currency="USD",
            steps=best_route.steps,
            default_step_yield=best_route.step_yield_assumed,
            target_scale_g=target_scale_g,
        )
        eff_moq_g = breakdown.moq_pack_g if breakdown.moq_pack_g is not None else record.moq_pack_g
        outlay_usd = (
            breakdown.min_outlay_cost.value
            if (breakdown.min_outlay_cost is not None and breakdown.min_outlay_complete)
            else (
                best_route.min_outlay_cost.value
                if (best_route.min_outlay_cost is not None and best_route.est_cost_complete)
                else None
            )
        )
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
            moq_pack_g=eff_moq_g,
            selected_pack_size=None,
            min_order_outlay_usd=outlay_usd,
            precursor_smiles=precursors,
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
        moq_pack_g=None,
        selected_pack_size=None,
        min_order_outlay_usd=None,
        precursor_smiles=[],
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
    "build_precursor_attribution_basis",
    "select_record_license",
    "to_procurement_constraint",
]
