"""Auto-Lab 15-repository wire node (`autolab-wire/v1`) for `buy-or-make` (Layer L5).

Deterministic offline make-vs-buy procurement triage enforcing `BuyOrMakeLean`
strict-null budget gating (`unpriced_never_passes_budget`), acyclic / unimolecular-hydrolytic
heavy-atom conservation, fail-closed packet/price input validation, explicit `data_source`
provenance (preventing synthetic reference stock from ever being confused with live vendor
catalogs), `labshift` (`L5`) release-time / MOQ scheduling handoff, and `lab-attrib` (`L4`)
precursor attribution basis export without live network requests or mutating working directories.
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
MAX_ALLOWED_WIRE_STEPS = 6

REPO_DIR = Path(__file__).resolve().parent
if str(REPO_DIR) not in sys.path:
    sys.path.insert(0, str(REPO_DIR))

from buyormake.handoff import (  # noqa: E402
    build_precursor_attribution_basis,
    to_procurement_constraint,
)
from buyormake.rank import evaluate_candidates, filter_decisions_from_ranked  # noqa: E402
from buyormake.retro_star import (  # noqa: E402
    RETRO_STAR_TEMPLATE_DELTAS,
    UNIMOLECULAR_RETRO_STAR_TEMPLATES,
)
from buyormake.schema import Route  # noqa: E402

# Deterministic offline reference catalog (strict `listed_catalog` USD/g offers, CC0-1.0,
# explicitly stamped `synthetic_reference_catalog` / `synthetic_demo`).
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

    Accepts bimolecular acyclic (`!@`) templates and registered unimolecular hydrolytic / FGI
    (`lactone_ring_opening`, `lactam_ring_opening`, etc. with `len(reactants) == 1`) templates,
    enforcing exact heavy-atom conservation across all chemical elements.
    """
    if not (1 <= route.n_steps <= max_steps):
        return False
    if not route.leaves or not all(route.leaves_purchasable):
        return False
    for step in route.steps:
        ref = step.template_or_reference or ""
        tpl_name = ref.split(":", 1)[0].strip()
        is_acyclic_bimol = "!@" in ref
        is_unimol_fgi = (
            tpl_name in UNIMOLECULAR_RETRO_STAR_TEMPLATES and len(step.reactants) == 1
        )
        if not (is_acyclic_bimol or is_unimol_fgi):
            return False
        delta = RETRO_STAR_TEMPLATE_DELTAS.get(tpl_name)
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


def _is_valid_finite_nonneg_number(val: Any) -> bool:
    return (
        not isinstance(val, bool)
        and isinstance(val, (int, float))
        and math.isfinite(float(val))
        and float(val) >= 0.0
    )


def _fail_closed_input_block(
    out: dict[str, Any],
    *,
    reason_detail: str,
    smiles_repr: str = "",
) -> dict[str, Any]:
    """Fail-closed wire response (`L5_PROCUREMENT_INVALID_INPUT`) without raising raw exceptions."""
    out["halted"] = True
    out["halt_reason"] = "L5_PROCUREMENT_INVALID_INPUT"
    out["procurement"] = {
        "smiles": smiles_repr,
        "status": "invalid",
        "listed_price_usd_per_g": None,
        "n_steps": None,
        "price_basis": "unknown",
        "within_budget": False,
        "data_source": "synthetic_reference_catalog",
        "catalog_origin": "synthetic_demo",
        "lead_time_days": None,
        "lead_time_complete": False,
        "release_time_minutes": None,
        "moq_pack_g": None,
        "selected_pack_size": None,
        "min_order_outlay_usd": None,
        "precursor_smiles": [],
        "precursor_basis": [],
        "reaction_steps": [],
        "validation_error": reason_detail,
    }
    out["trace"].append(
        {
            "repo": REPO_ID,
            "layer": LAYER,
            "status": "BLOCK",
            "summary": f"L5 input validation blocked packet: {reason_detail}",
        }
    )
    out["bottom_reasoning_traces"].append(
        {
            "repo": REPO_ID,
            "layer": LAYER,
            "invariant": "unpriced_never_passes_budget",
            "first_principles": (
                "Malformed, non-finite, negative, boolean, or unparseable procurement inputs "
                "fail closed at Layer L5 and never pass budget or substitute default molecules."
            ),
        }
    )
    return out


def process_packet(packet: dict[str, Any]) -> dict[str, Any]:
    """Evaluate make-vs-buy procurement triage on an `autolab-wire/v1` packet."""
    if not isinstance(packet, dict):
        fresh: dict[str, Any] = {
            "schema_version": WIRE_SCHEMA_VERSION,
            "halted": False,
            "halt_reason": None,
            "trace": [],
            "bottom_reasoning_traces": [],
        }
        return _fail_closed_input_block(
            fresh, reason_detail=f"packet must be a dict, got {type(packet).__name__}"
        )

    out: dict[str, Any] = copy.deepcopy(packet)
    out.setdefault("schema_version", WIRE_SCHEMA_VERSION)
    if "halted" not in out:
        out["halted"] = False
    out.setdefault("halt_reason", None)

    if not isinstance(out.get("trace", []), list) or not isinstance(
        out.get("bottom_reasoning_traces", []), list
    ):
        out["trace"] = out["trace"] if isinstance(out.get("trace"), list) else []
        out["bottom_reasoning_traces"] = (
            out["bottom_reasoning_traces"]
            if isinstance(out.get("bottom_reasoning_traces"), list)
            else []
        )
        return _fail_closed_input_block(
            out, reason_detail="packet trace / bottom_reasoning_traces must be lists"
        )
    out.setdefault("trace", [])
    out.setdefault("bottom_reasoning_traces", [])

    if not isinstance(out["halted"], bool):
        return _fail_closed_input_block(
            out, reason_detail=f"packet halted must be a bool, got {type(out['halted']).__name__}"
        )

    if out["halted"] is True:
        return out

    # Strict candidate_smiles extraction (never silently fall back to aspirin when explicitly malformed!)
    if "candidate_smiles" in out:
        raw_smi = out["candidate_smiles"]
        if not isinstance(raw_smi, str) or not raw_smi.strip():
            return _fail_closed_input_block(
                out,
                reason_detail="candidate_smiles must be a non-empty SMILES string",
                smiles_repr="" if raw_smi is None else str(raw_smi),
            )
        smiles = raw_smi.strip()
    elif "decision" in out:
        dec_obj = out["decision"]
        if dec_obj is not None and not isinstance(dec_obj, dict):
            return _fail_closed_input_block(
                out, reason_detail="decision field must be a dict when provided"
            )
        if isinstance(dec_obj, dict) and "candidate_smiles" in dec_obj:
            raw_smi = dec_obj["candidate_smiles"]
            if not isinstance(raw_smi, str) or not raw_smi.strip():
                return _fail_closed_input_block(
                    out,
                    reason_detail="decision.candidate_smiles must be a non-empty SMILES string",
                    smiles_repr="" if raw_smi is None else str(raw_smi),
                )
            smiles = raw_smi.strip()
        else:
            smiles = DEFAULT_CANDIDATE_SMILES
    else:
        smiles = DEFAULT_CANDIDATE_SMILES

    # Validate budget_usd_per_g
    if "budget_usd_per_g" in out:
        raw_budget = out["budget_usd_per_g"]
        if not _is_valid_finite_nonneg_number(raw_budget):
            return _fail_closed_input_block(
                out,
                reason_detail=f"budget_usd_per_g must be a finite number >= 0, got {raw_budget!r}",
                smiles_repr=smiles,
            )
        budget = float(raw_budget)
    else:
        budget = DEFAULT_BUDGET_USD_PER_G

    # Validate max_steps
    if "max_steps" in out:
        raw_ms = out["max_steps"]
        if (
            isinstance(raw_ms, bool)
            or not isinstance(raw_ms, int)
            or not (0 <= raw_ms <= MAX_ALLOWED_WIRE_STEPS)
        ):
            return _fail_closed_input_block(
                out,
                reason_detail=(
                    f"max_steps must be an integer in [0, {MAX_ALLOWED_WIRE_STEPS}], got {raw_ms!r}"
                ),
                smiles_repr=smiles,
            )
        max_steps = int(raw_ms)
    else:
        max_steps = 2

    # Validate require_priced_within_budget
    if "require_priced_within_budget" in out:
        raw_req = out["require_priced_within_budget"]
        if not isinstance(raw_req, bool):
            return _fail_closed_input_block(
                out,
                reason_detail=(
                    f"require_priced_within_budget must be a bool, got {raw_req!r}"
                ),
                smiles_repr=smiles,
            )
        require_priced = raw_req
    else:
        require_priced = False

    # Validate optional planner, step_yield, target_scale_g, max_outlay_usd, and upstream price hints
    planner = out.get("planner", "lookup_smarts")
    if not isinstance(planner, str) or planner not in {"lookup_smarts", "retro_star"}:
        return _fail_closed_input_block(
            out,
            reason_detail=f"planner must be 'lookup_smarts' or 'retro_star', got {planner!r}",
            smiles_repr=smiles,
        )

    if "step_yield" in out:
        raw_yield = out["step_yield"]
        if (
            isinstance(raw_yield, bool)
            or not isinstance(raw_yield, (int, float))
            or not math.isfinite(float(raw_yield))
            or not (0.0 < float(raw_yield) <= 1.0)
        ):
            return _fail_closed_input_block(
                out,
                reason_detail=f"step_yield must be a finite float in (0, 1], got {raw_yield!r}",
                smiles_repr=smiles,
            )
        step_yield = float(raw_yield)
    else:
        step_yield = 1.0

    if "target_scale_g" in out:
        raw_scale = out["target_scale_g"]
        if (
            isinstance(raw_scale, bool)
            or not isinstance(raw_scale, (int, float))
            or not math.isfinite(float(raw_scale))
            or float(raw_scale) <= 0.0
        ):
            return _fail_closed_input_block(
                out,
                reason_detail=f"target_scale_g must be a finite float > 0, got {raw_scale!r}",
                smiles_repr=smiles,
            )
        target_scale_g = float(raw_scale)
    else:
        target_scale_g = 1.0

    max_outlay_usd: float | None = None
    if "max_outlay_usd" in out and out["max_outlay_usd"] is not None:
        raw_outlay = out["max_outlay_usd"]
        if not _is_valid_finite_nonneg_number(raw_outlay):
            return _fail_closed_input_block(
                out,
                reason_detail=f"max_outlay_usd must be a finite number >= 0, got {raw_outlay!r}",
                smiles_repr=smiles,
            )
        max_outlay_usd = float(raw_outlay)

    if "listed_price_usd_per_g" in out and out["listed_price_usd_per_g"] is not None:
        if not _is_valid_finite_nonneg_number(out["listed_price_usd_per_g"]):
            return _fail_closed_input_block(
                out,
                reason_detail=(
                    f"listed_price_usd_per_g must be a finite number >= 0, "
                    f"got {out['listed_price_usd_per_g']!r}"
                ),
                smiles_repr=smiles,
            )

    custom_stock_paths: list[Path] | None = None
    if "stock_paths" in out and out["stock_paths"] is not None:
        raw_sps = out["stock_paths"]
        if not isinstance(raw_sps, list) or not raw_sps or not all(
            isinstance(p, (str, Path)) and Path(p).exists() for p in raw_sps
        ):
            return _fail_closed_input_block(
                out,
                reason_detail="stock_paths must be a non-empty list of existing file paths",
                smiles_repr=smiles,
            )
        custom_stock_paths = [Path(p) for p in raw_sps]

    with tempfile.TemporaryDirectory(prefix="buyormake_wire_") as tmp_raw:
        tmp_dir = Path(tmp_raw)
        if custom_stock_paths is not None:
            active_stock_paths = custom_stock_paths
            using_synthetic_ref = False
        else:
            stock_path = tmp_dir / "reference_stock.tsv"
            stock_path.write_text(_REFERENCE_STOCK_TSV, encoding="utf-8")
            active_stock_paths = [stock_path]
            using_synthetic_ref = True

        ranked, records = evaluate_candidates(
            [smiles],
            max_steps=max_steps,
            backends=("local_stock",),
            route_backend="template",
            planner=planner,
            offline=True,
            cache_path=tmp_dir / "wire_cache.sqlite",
            stock_paths=active_stock_paths,
            reaction_db_path=tmp_dir / "no_reactions.sqlite",
            scscore_weights=REPO_DIR / "data" / "scscore" / "model.ckpt-10654.as_numpy.json.gz",
            include_scores=True,
            desalt_fallback=True,
            default_step_yield=step_yield,
            target_scale_g=target_scale_g,
            currency="USD",
        )
        filter_summary = filter_decisions_from_ranked(
            ranked,
            budget=budget,
            budget_currency="USD",
            require_complete_cost=True,
            max_outlay_budget=max_outlay_usd,
        )

    res = ranked[0]
    if res.status == "invalid":
        return _fail_closed_input_block(
            out,
            reason_detail=f"unparseable SMILES {smiles!r}",
            smiles_repr=smiles,
        )

    decision = filter_summary.decisions[0]
    rec = records[0] if records else None
    constraint = (
        to_procurement_constraint(rec, target_scale_g=target_scale_g)
        if rec is not None
        else None
    )

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
        and decision.keep
    )
    if max_outlay_usd is not None and (
        min_order_outlay_usd is None or min_order_outlay_usd > max_outlay_usd
    ):
        within_budget = False

    if require_priced and (price is None or not within_budget or price > budget):
        out["halted"] = True
        out["halt_reason"] = "L5_PROCUREMENT_BUDGET_BLOCK"

    data_source = (
        "synthetic_reference_catalog"
        if using_synthetic_ref
        else str(res.data_source)
    )
    catalog_origin = (
        "synthetic_demo"
        if using_synthetic_ref
        else (str(rec.origin) if rec is not None else "real")
    )

    out["procurement"] = {
        "smiles": smiles,
        "status": status,
        "listed_price_usd_per_g": price,
        "n_steps": n_steps,
        "price_basis": price_basis,
        "within_budget": within_budget,
        "data_source": data_source,
        "catalog_origin": catalog_origin,
        "planner": planner,
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
                f"Triage {smiles}: status={status}, price={price_str} ({price_basis}, "
                f"source={data_source}), n_steps={n_steps}, budget=${budget:.2f}/g, "
                f"within_budget={within_budget}{release_str}"
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
    "MAX_ALLOWED_WIRE_STEPS",
    "REPO_DIR",
    "REPO_ID",
    "WIRE_SCHEMA_VERSION",
    "process_packet",
]
