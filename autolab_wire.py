"""Auto-Lab 15-repository wire node (`autolab-wire/v1`) for `buy-or-make` (Layer L5).

Deterministic offline make-vs-buy procurement triage enforcing `BuyOrMakeLean`
strict-null budget gating (`unpriced_never_passes_budget`), acyclic / unimolecular-hydrolytic
heavy-atom conservation, `labshift` (`L5`) release-time / MOQ scheduling handoff, and
`lab-attrib` (`L4`) precursor attribution basis export without live network requests or
mutating working directories.
"""

from __future__ import annotations

import copy
import math
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

from rdkit import Chem

REPO_ID = "buy-or-make"
LAYER = "L5"
WIRE_SCHEMA_VERSION = "autolab-wire/v1"
DEFAULT_CANDIDATE_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
DEFAULT_BUDGET_USD_PER_G = 25.0

REPO_DIR = Path(__file__).resolve().parent
if str(REPO_DIR) not in sys.path:
    sys.path.insert(0, str(REPO_DIR))

from buyormake.handoff import (  # noqa: E402
    build_precursor_attribution_basis,
    to_procurement_constraint,
)
from buyormake.rank import evaluate_candidates, filter_decisions_from_ranked  # noqa: E402
from buyormake.routes import (  # noqa: E402
    RETRO_TEMPLATE_DELTAS,
    UNIMOLECULAR_RING_OPENING_TEMPLATES,
)
from buyormake.schema import Route  # noqa: E402

# Deterministic offline reference catalog (strict `listed_catalog` USD/g offers, CC0-1.0).
_REFERENCE_STOCK_ROWS: tuple[str, ...] = (
    "smiles\tvendor\tcatalog_id\tpack_value\tpack_unit\tprice\tcurrency\tlead_time_days\tin_stock\tprice_basis\tsource_id\tsource_url\tlicense",
    "CC(=O)Oc1ccccc1C(=O)O\tautolab_ref_catalog\tREF-ASP-001\t1.0\tg\t0.42\tUSD\t2\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-ASP-001\tCC0-1.0",
    "CC(=O)Nc1ccc(O)cc1\tautolab_ref_catalog\tREF-PAR-002\t1.0\tg\t0.38\tUSD\t2\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-PAR-002\tCC0-1.0",
    "CC(C)Cc1ccc(C(C)C(=O)O)cc1\tautolab_ref_catalog\tREF-IBU-003\t1.0\tg\t0.65\tUSD\t2\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-IBU-003\tCC0-1.0",
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O\tautolab_ref_catalog\tREF-CAF-004\t1.0\tg\t0.55\tUSD\t2\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-CAF-004\tCC0-1.0",
    "O=C(O)c1ccccc1O\tautolab_ref_catalog\tREF-SAL-005\t1.0\tg\t0.18\tUSD\t1\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-SAL-005\tCC0-1.0",
    "CC(=O)O\tautolab_ref_catalog\tREF-ACE-006\t1.0\tg\t0.08\tUSD\t1\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-ACE-006\tCC0-1.0",
    "O=C(O)c1ccccc1\tautolab_ref_catalog\tREF-BZA-007\t1.0\tg\t0.15\tUSD\t1\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-BZA-007\tCC0-1.0",
    "Nc1ccccc1\tautolab_ref_catalog\tREF-ANI-008\t1.0\tg\t0.20\tUSD\t1\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-ANI-008\tCC0-1.0",
    "COc1ccc(C(=O)O)cc1\tautolab_ref_catalog\tREF-MBA-009\t1.0\tg\t0.40\tUSD\t2\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-MBA-009\tCC0-1.0",
    "Nc1ccc(Br)cc1\tautolab_ref_catalog\tREF-BRA-010\t1.0\tg\t0.80\tUSD\t2\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-BRA-010\tCC0-1.0",
    "OB(O)c1ccccn1\tautolab_ref_catalog\tREF-PYB-011\t1.0\tg\t3.50\tUSD\t3\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-PYB-011\tCC0-1.0",
    "O=S(=O)(Cl)c1ccc(F)cc1\tautolab_ref_catalog\tREF-FBS-012\t1.0\tg\t1.20\tUSD\t2\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-FBS-012\tCC0-1.0",
    "CCO\tautolab_ref_catalog\tREF-ETH-013\t1.0\tg\t0.05\tUSD\t1\ttrue\tlisted_catalog\tlocal_stock\thttps://autolab.local/catalog/REF-ETH-013\tCC0-1.0",
)
_REFERENCE_STOCK_TSV = "\n".join(_REFERENCE_STOCK_ROWS) + "\n"


def _heavy_counts(smiles: str) -> Counter[str] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def _route_is_acyclic_and_atom_balanced(route: Route, max_steps: int) -> bool:
    """Verify `BuyOrMakeLean.routeStoichiometryOK` on every retrosynthetic step.

    Accepts bimolecular acyclic (`!@`) templates and registered unimolecular hydrolytic
    ring-opening (`lactone_ring_opening`, `lactam_ring_opening` with `len(reactants) == 1`)
    templates, enforcing exact heavy-atom conservation across all chemical elements.
    """
    if not (1 <= route.n_steps <= max_steps):
        return False
    if not route.leaves or not all(route.leaves_purchasable):
        return False
    for step in route.steps:
        ref = step.template_or_reference or ""
        tpl_name = ref.split(":", 1)[0].strip()
        is_acyclic_bimol = "!@" in ref
        is_unimol_ring_open = (
            tpl_name in UNIMOLECULAR_RING_OPENING_TEMPLATES and len(step.reactants) == 1
        )
        if not (is_acyclic_bimol or is_unimol_ring_open):
            return False
        delta = RETRO_TEMPLATE_DELTAS.get(tpl_name)
        prod_counts = _heavy_counts(step.product)
        if delta is None or prod_counts is None or not step.reactants:
            return False
        expected = Counter(prod_counts)
        for sym, inc in delta.items():
            expected[sym] += inc
            if expected[sym] <= 0:
                del expected[sym]
        reactant_counts: Counter[str] = Counter()
        for r_smi in step.reactants:
            rc = _heavy_counts(r_smi)
            if rc is None:
                return False
            reactant_counts.update(rc)
        if reactant_counts != expected:
            return False
    return True


def process_packet(packet: dict[str, Any]) -> dict[str, Any]:
    """Evaluate make-vs-buy procurement triage on an `autolab-wire/v1` packet."""
    out: dict[str, Any] = copy.deepcopy(packet)
    out.setdefault("schema_version", WIRE_SCHEMA_VERSION)
    out.setdefault("halted", False)
    out.setdefault("halt_reason", None)
    out.setdefault("trace", [])
    out.setdefault("bottom_reasoning_traces", [])

    if out["halted"] is True:
        return out

    if "candidate_smiles" in out and out["candidate_smiles"] is not None:
        smiles = str(out["candidate_smiles"])
    elif isinstance(out.get("decision"), dict) and out["decision"].get("candidate_smiles") is not None:
        smiles = str(out["decision"]["candidate_smiles"])
    else:
        smiles = DEFAULT_CANDIDATE_SMILES

    budget = float(out.get("budget_usd_per_g", DEFAULT_BUDGET_USD_PER_G))
    max_steps = int(out.get("max_steps", 2))
    require_priced = bool(out.get("require_priced_within_budget", False))

    with tempfile.TemporaryDirectory(prefix="buyormake_wire_") as tmp_raw:
        tmp_dir = Path(tmp_raw)
        stock_path = tmp_dir / "reference_stock.tsv"
        stock_path.write_text(_REFERENCE_STOCK_TSV, encoding="utf-8")

        ranked, records = evaluate_candidates(
            [smiles],
            max_steps=max_steps,
            backends=("local_stock",),
            route_backend="template",
            offline=True,
            cache_path=tmp_dir / "wire_cache.sqlite",
            stock_paths=[stock_path],
            reaction_db_path=tmp_dir / "no_reactions.sqlite",
            scscore_weights=REPO_DIR / "data" / "scscore" / "model.ckpt-10654.as_numpy.json.gz",
            include_scores=True,
            desalt_fallback=True,
            currency="USD",
        )
        filter_summary = filter_decisions_from_ranked(
            ranked,
            budget=budget if math.isfinite(budget) else None,
            budget_currency="USD",
            require_complete_cost=True,
        )

    res = ranked[0]
    decision = filter_summary.decisions[0]
    rec = records[0] if records else None
    constraint = to_procurement_constraint(rec) if rec is not None else None

    status = res.status
    n_steps = res.n_steps
    price: float | None = None
    price_basis = "unknown"
    lead_time_days: int | None = constraint.lead_time_days if constraint is not None else None
    lead_time_complete: bool = constraint.lead_time_complete if constraint is not None else False
    release_time_minutes: int | None = (
        constraint.release_time_minutes if constraint is not None else None
    )
    moq_pack_g: float | None = constraint.moq_pack_g if constraint is not None else None
    selected_pack_size: dict[str, Any] | None = (
        constraint.selected_pack_size.model_dump(mode="json", exclude_none=True)
        if (constraint is not None and constraint.selected_pack_size is not None)
        else None
    )
    min_order_outlay_usd: float | None = (
        constraint.min_order_outlay_usd if constraint is not None else None
    )
    precursor_smiles: list[str] = (
        list(constraint.precursor_smiles) if constraint is not None else []
    )
    precursor_basis: list[dict[str, Any]] = []
    reaction_steps: list[dict[str, Any]] = []

    if status == "buyable":
        if (
            res.best_price is not None
            and res.best_price.value is not None
            and res.best_price.unit == "USD/g"
        ):
            price = round(float(res.best_price.value), 4)
            price_basis = "listed_catalog"
    elif status == "makeable":
        solved_routes = (
            [r for r in rec.routes if r.leaves and all(r.leaves_purchasable)]
            if rec is not None
            else []
        )
        if not solved_routes or not _route_is_acyclic_and_atom_balanced(
            solved_routes[0], max_steps=max_steps
        ):
            status = "unknown"
            n_steps = None
            lead_time_days = None
            lead_time_complete = False
            release_time_minutes = None
            moq_pack_g = None
            selected_pack_size = None
            min_order_outlay_usd = None
            precursor_smiles = []
        else:
            best_route = solved_routes[0]
            if rec is not None:
                precursor_basis, reaction_steps = build_precursor_attribution_basis(
                    best_route, rec.leaf_offers, currency="USD"
                )
            if (
                res.cost_complete
                and res.est_cost is not None
                and res.est_cost.value is not None
                and res.est_cost.unit in {"USD", "USD/g"}
            ):
                price = round(float(res.est_cost.value), 4)
                price_basis = "listed_catalog"

    within_budget = bool(
        price is not None
        and math.isfinite(price)
        and price <= budget
        and (decision.keep if math.isfinite(budget) else True)
    )

    if require_priced and (price is None or not within_budget or price > budget):
        out["halted"] = True
        out["halt_reason"] = "L5_PROCUREMENT_BUDGET_BLOCK"

    out["procurement"] = {
        "smiles": smiles,
        "status": status,
        "listed_price_usd_per_g": price,
        "n_steps": n_steps,
        "price_basis": price_basis,
        "within_budget": within_budget,
        "lead_time_days": lead_time_days,
        "lead_time_complete": lead_time_complete,
        "release_time_minutes": release_time_minutes,
        "moq_pack_g": moq_pack_g,
        "selected_pack_size": selected_pack_size,
        "min_order_outlay_usd": min_order_outlay_usd,
        "precursor_smiles": precursor_smiles,
        "precursor_basis": precursor_basis,
        "reaction_steps": reaction_steps,
    }

    price_str = f"${price:.2f}/g" if price is not None else "null"
    release_str = (
        f", release_min={release_time_minutes}"
        if release_time_minutes is not None
        else ""
    )
    out["trace"].append(
        {
            "repo": REPO_ID,
            "layer": LAYER,
            "status": "BLOCK" if out["halted"] else "PASS",
            "summary": (
                f"Triage {smiles}: status={status}, price={price_str} ({price_basis}), "
                f"n_steps={n_steps}, budget=${budget:.2f}/g, within_budget={within_budget}{release_str}"
            ),
        }
    )
    out["bottom_reasoning_traces"].append(
        {
            "repo": REPO_ID,
            "layer": LAYER,
            "invariant": "unpriced_never_passes_budget",
            "first_principles": (
                "Missing or incomplete catalog prices stay null and never pass a numeric USD/g "
                "budget gate; synthesis routes must conserve heavy-atom counts across every "
                "element over purchasable leaves without atom-duplicating ring cleavage."
            ),
        }
    )
    return out


__all__ = [
    "DEFAULT_BUDGET_USD_PER_G",
    "DEFAULT_CANDIDATE_SMILES",
    "LAYER",
    "REPO_DIR",
    "REPO_ID",
    "WIRE_SCHEMA_VERSION",
    "process_packet",
]
