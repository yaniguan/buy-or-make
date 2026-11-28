"""Unified schema for buy-or-make (Pydantic v2, extra='forbid').

Field names follow the assayfold hub schema (Source, Quantity, Uncertainty,
AgentRef, SoftwareRef, Provenance, Molecule, origin, schema_version, record_hash)
so availability records can be exported to an open-data hub without renaming.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Literal, Optional

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .canonical import canonical_json, sha256_hex

SCHEMA_VERSION = "0.1.0"

# Conservative license restrictiveness ranking (higher = more restrictive; wins when combining).
_LICENSE_RESTRICTIVENESS: dict[str, int] = {
    "CC0-1.0": 10,
    "MIT": 20,
    "BSD-2-Clause": 30,
    "BSD-3-Clause": 30,
    "Apache-2.0": 35,
    "CC-BY-3.0": 50,
    "CC-BY-4.0": 50,
    "CC-BY-SA-3.0": 70,
    "CC-BY-SA-4.0": 70,
    "CC-BY-NC-4.0": 85,
    "CC-BY-NC-SA-4.0": 90,
    "LicenseRef-Custom-Terms": 100,
    "custom/terms": 100,
}


def select_record_license(licenses: Iterable[str], default: str = "MIT") -> str:
    """Select the single governing license for a combined record (most restrictive wins)."""
    cleaned = sorted({lic.strip() for lic in licenses if lic and lic.strip()})
    if not cleaned:
        return default
    return max(cleaned, key=lambda lic: (_LICENSE_RESTRICTIVENESS.get(lic, 95), lic))


def make_assayfold_record_id(source_id: str, dataset_id: str, native_id: str) -> str:
    """Deterministic assayfold record_id ('af:' + 24 hex chars of SHA-256)."""
    return "af:" + sha256_hex(f"{source_id}|{dataset_id}|{native_id}".encode("utf-8"))[:24]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid")


# ---------- shared value types (field names follow the assayfold hub schema) ----------


class Uncertainty(_Model):
    kind: Literal["std", "interval", "none"] = "none"
    std: Optional[float] = None
    low: Optional[float] = None
    high: Optional[float] = None
    n: Optional[int] = None
    method: Optional[str] = None

    @model_validator(mode="after")
    def _validate_uncertainty(self) -> "Uncertainty":
        for fname in ("std", "low", "high"):
            val = getattr(self, fname)
            if val is not None and not math.isfinite(val):
                raise ValueError(f"Uncertainty.{fname} must be finite")
        if self.std is not None and self.std < 0:
            raise ValueError("Uncertainty.std must be >= 0")
        if self.n is not None and self.n < 0:
            raise ValueError("Uncertainty.n must be >= 0")
        if self.low is not None and self.high is not None and self.low > self.high:
            raise ValueError("Uncertainty.low must be <= Uncertainty.high")
        return self


class Quantity(_Model):
    value: Optional[float] = None
    unit: str = "1"
    raw_unit_str: Optional[str] = None
    uncertainty: Optional[Uncertainty] = None

    @model_validator(mode="after")
    def _validate_quantity(self) -> "Quantity":
        if self.value is not None and not math.isfinite(self.value):
            raise ValueError("Quantity.value must be a finite float")
        if not self.unit or not self.unit.strip():
            raise ValueError("Quantity.unit must be non-empty")
        return self


class BlobRef(_Model):
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    media_type: str
    size_bytes: int = Field(ge=0)
    filename_orig: Optional[str] = None
    url: Optional[str] = None
    role: Literal["raw", "derived", "thumbnail"] = "raw"


class AgentRef(_Model):
    type: Literal["human", "robot", "model", "software", "instrument"] = "software"
    id: str
    name: Optional[str] = None
    version: Optional[str] = None


class SoftwareRef(_Model):
    name: str
    version: Optional[str] = None


class Provenance(_Model):
    agent: AgentRef = Field(default_factory=lambda: AgentRef(id="buy-or-make", name="buy-or-make"))
    recorded_utc: str = Field(default_factory=utc_now)
    source_url: Optional[str] = None
    source_file_sha256: Optional[str] = None
    software: list[SoftwareRef] = Field(default_factory=list)
    used: list[str] = Field(default_factory=list)
    derived_from: list[str] = Field(default_factory=list)


class Source(_Model):
    id: str = Field(pattern=r"^[a-z0-9_]+$")
    name: str
    homepage: str
    license: str
    spdx_license: Optional[str] = None
    license_url: Optional[str] = None
    license_evidence: str
    license_quote: Optional[str] = None
    license_checked: str = Field(pattern=r"^\d{4}-\d{2}-\d{2}$")
    license_check_method: Literal["direct", "indirect"] = "direct"
    commercial_use_ok: bool
    redistribution_ok: bool
    attribution: str
    modalities: list[str] = Field(default_factory=list)
    kinds: list[str] = Field(default_factory=list)
    access: Literal["api", "download", "huggingface", "github", "builtin"]
    size_full: Optional[str] = None
    status: Literal["integrated", "fetch_only", "not_used", "planned"] = "integrated"
    notes: Optional[str] = None

    @model_validator(mode="after")
    def _populate_spdx_license(self) -> "Source":
        if self.spdx_license is not None and not self.spdx_license.strip():
            raise ValueError("Source.spdx_license must be non-empty when provided")
        if self.spdx_license is None:
            self.spdx_license = (
                "LicenseRef-Custom-Terms" if self.license == "custom/terms" else self.license
            )
        return self


class Molecule(_Model):
    kind: Literal["molecule"] = "molecule"
    smiles: str
    inchikey: Optional[str] = None
    name: Optional[str] = None
    formula: Optional[str] = None


# ---------- buy-or-make domain types ----------


class Offer(_Model):
    inchikey: str
    vendor: str
    catalog_id: Optional[str] = None
    pack_size: Optional[Quantity] = None
    purity: Optional[str] = None
    price: Optional[float] = None
    currency: Optional[str] = None
    price_basis: Literal["listed_catalog", "unknown"] = "unknown"
    lead_time_days: Optional[int] = None
    in_stock: Optional[bool] = None
    source_id: str
    source_url: str
    fetched_utc: str = Field(default_factory=utc_now)
    license: str

    @model_validator(mode="after")
    def _validate_offer(self) -> "Offer":
        if self.price is not None:
            if isinstance(self.price, bool) or not math.isfinite(self.price) or self.price < 0:
                raise ValueError("Offer.price must be a finite number >= 0")
            if self.currency is None or not self.currency.strip():
                raise ValueError("Offer.currency is required whenever Offer.price is set")
        else:
            if self.currency is not None:
                raise ValueError("Offer.currency must be None when Offer.price is None")
            if self.price_basis == "listed_catalog":
                raise ValueError("Offer.price_basis='listed_catalog' requires a non-null price")
        if self.lead_time_days is not None:
            if isinstance(self.lead_time_days, bool) or self.lead_time_days < 0:
                raise ValueError("Offer.lead_time_days must be >= 0")
        if self.pack_size is not None and self.pack_size.value is not None:
            if self.pack_size.value <= 0:
                raise ValueError("Offer.pack_size.value must be > 0")
        return self


class RouteStep(_Model):
    product: str
    reactants: list[str]
    reaction_smiles: Optional[str] = None
    template_or_reference: str
    source: str


class Route(_Model):
    target_inchikey: str
    n_steps: int
    tree_depth: Optional[int] = None
    steps: list[RouteStep]
    leaves: list[str]
    leaves_purchasable: list[bool]
    method: str
    score: Optional[float] = None
    est_cost: Optional[Quantity] = Field(
        default=None,
        description=(
            "Estimated stoichiometric leaf-only cost per 1 g of target in '<currency>/g' "
            "(e.g. 'USD/g'), independent of target_scale_g."
        ),
    )
    est_cost_method: str
    est_cost_complete: bool
    est_cost_yield_adjusted: Optional[Quantity] = Field(
        default=None,
        description=(
            "Step-depth yield-compounded leaf-only cost per 1 g of target in '<currency>/g' "
            "(e.g. 'USD/g')."
        ),
    )
    min_outlay_cost: Optional[Quantity] = Field(
        default=None,
        description=(
            "Minimum discrete pack MOQ cash outlay for the full target_scale_g order in "
            "absolute '<currency>' (e.g. 'USD')."
        ),
    )
    step_yield_assumed: float = 1.0

    @model_validator(mode="after")
    def _validate_route(self) -> "Route":
        if self.n_steps < 0:
            raise ValueError("Route.n_steps must be >= 0")
        if self.tree_depth is not None and (isinstance(self.tree_depth, bool) or self.tree_depth < 0):
            raise ValueError("Route.tree_depth must be >= 0 when provided")
        if len(self.leaves) != len(self.leaves_purchasable):
            raise ValueError("Route.leaves and Route.leaves_purchasable must have equal length")
        if self.score is not None and not math.isfinite(self.score):
            raise ValueError("Route.score must be finite")
        if (
            isinstance(self.step_yield_assumed, bool)
            or not math.isfinite(self.step_yield_assumed)
            or not (0.0 < self.step_yield_assumed <= 1.0)
        ):
            raise ValueError("Route.step_yield_assumed must be in (0, 1]")
        return self


class SAScores(_Model):
    sascore: Optional[float] = None
    scscore: Optional[float] = None
    rascore: Optional[float] = None

    @model_validator(mode="after")
    def _validate_scores(self) -> "SAScores":
        for fname in ("sascore", "scscore", "rascore"):
            val = getattr(self, fname)
            if val is not None and not math.isfinite(val):
                raise ValueError(f"SAScores.{fname} must be finite")
        return self


OriginLiteral = Literal["real", "simulated", "model_generated", "synthetic_demo"]


class AssayfoldRecord(_Model):
    """Top-level `assayfold` hub `Record(kind='molecule')` envelope."""

    record_id: str
    source_id: str = Field(pattern=r"^[a-z0-9_]+$")
    dataset_id: str
    native_id: str
    kind: Literal["molecule"] = "molecule"
    parent_id: Optional[str] = None
    text: Optional[str] = None
    payload: Molecule
    properties: dict[str, Quantity] = Field(default_factory=dict)
    blobs: list[BlobRef] = Field(default_factory=list)
    license: str
    origin: OriginLiteral = "real"
    provenance: Provenance = Field(default_factory=Provenance)
    tags: list[str] = Field(default_factory=list)
    schema_version: str = SCHEMA_VERSION
    record_hash: Optional[str] = None

    @model_validator(mode="after")
    def _kind_matches(self) -> "AssayfoldRecord":
        if self.kind != self.payload.kind:
            raise ValueError(f"kind {self.kind!r} != payload kind {self.payload.kind!r}")
        return self

    def content_hash(self) -> str:
        """SHA-256 over canonical JSON. Matches `assayfold.schema.Record.content_hash()`."""
        d = self.model_dump(mode="json", exclude={"record_hash"})
        d["provenance"].pop("recorded_utc", None)
        return sha256_hex(canonical_json(d))

    def sealed(self) -> "AssayfoldRecord":
        return self.model_copy(update={"record_hash": self.content_hash()})


InChIKeyMatchLiteral = Literal["exact", "salt_equivalent", "stereoisomer_relaxed", "none"]
DataSourceLiteral = Literal[
    "live_vendor_registry",
    "local_vendor_catalog",
    "open_stopping_stock",
    "synthetic_reference_catalog",
    "unmatched",
    "unknown",
]


class AvailabilityRecord(_Model):
    molecule: Molecule
    offers: list[Offer] = Field(default_factory=list)
    routes: list[Route] = Field(default_factory=list)
    leaf_offers: dict[str, list[Offer]] = Field(default_factory=dict)
    lead_time_days: Optional[int] = None
    moq_pack_g: Optional[float] = None
    precursor_smiles: list[str] = Field(default_factory=list)
    inchikey_match: InChIKeyMatchLiteral = "none"
    matched_offer_inchikeys: list[str] = Field(default_factory=list)
    stereoisomer_offers: list[Offer] = Field(default_factory=list)
    sa_scores: SAScores = Field(default_factory=SAScores)
    data_source: DataSourceLiteral = "unknown"
    origin: OriginLiteral = "real"
    provenance: Provenance = Field(default_factory=Provenance)
    licenses: list[str] = Field(default_factory=list)
    schema_version: str = SCHEMA_VERSION
    record_hash: Optional[str] = None

    @model_validator(mode="after")
    def _validate_availability_record(self) -> "AvailabilityRecord":
        if self.lead_time_days is not None:
            if isinstance(self.lead_time_days, bool) or self.lead_time_days < 0:
                raise ValueError("AvailabilityRecord.lead_time_days must be >= 0")
        if self.moq_pack_g is not None:
            if (
                isinstance(self.moq_pack_g, bool)
                or not math.isfinite(self.moq_pack_g)
                or self.moq_pack_g <= 0
            ):
                raise ValueError("AvailabilityRecord.moq_pack_g must be a finite number > 0")
        return self

    def content_hash(self) -> str:
        """SHA-256 over canonical JSON. Excludes record_hash and provenance.recorded_utc."""
        d = self.model_dump(mode="json", exclude={"record_hash"})
        d["provenance"].pop("recorded_utc", None)
        return sha256_hex(canonical_json(d))

    def sealed(self) -> "AvailabilityRecord":
        """Return a copy with record_hash filled in."""
        return self.model_copy(update={"record_hash": self.content_hash()})

    def to_assayfold_record(
        self,
        *,
        source_id: str = "buy_or_make",
        dataset_id: str = "availability_triage",
        native_id: Optional[str] = None,
    ) -> AssayfoldRecord:
        """Convert this `AvailabilityRecord` into a sealed `assayfold` `Record(kind='molecule')`."""
        from .cost import cheapest_offer_for_target, purchasable_offers

        nid = str(native_id or self.molecule.inchikey or self.molecule.smiles)
        purchasable = purchasable_offers(self.offers)
        props: dict[str, Quantity] = {
            "n_offers": Quantity(value=float(len(self.offers)), unit="1"),
        }
        best_price = cheapest_offer_for_target(purchasable if purchasable else self.offers)
        if best_price is not None:
            props["best_price"] = best_price
        if self.lead_time_days is not None:
            props["lead_time_days"] = Quantity(value=float(self.lead_time_days), unit="d")
        if self.moq_pack_g is not None:
            props["moq_pack_g"] = Quantity(value=float(self.moq_pack_g), unit="g")

        solved = [r for r in self.routes if r.leaves and all(r.leaves_purchasable)]
        if purchasable:
            status_tag = "status:buyable"
            props["min_synthesis_steps"] = Quantity(value=0.0, unit="1")
        elif solved:
            status_tag = "status:makeable"
            best_route = solved[0]
            props["min_synthesis_steps"] = Quantity(value=float(best_route.n_steps), unit="1")
            if best_route.est_cost is not None:
                props["est_route_cost"] = best_route.est_cost
            if best_route.score is not None:
                props["route_score"] = Quantity(value=float(best_route.score), unit="1")
        else:
            status_tag = "status:unknown"

        if self.sa_scores.sascore is not None:
            props["sascore"] = Quantity(value=self.sa_scores.sascore, unit="1")
        if self.sa_scores.scscore is not None:
            props["scscore"] = Quantity(value=self.sa_scores.scscore, unit="1")
        if self.sa_scores.rascore is not None:
            props["rascore"] = Quantity(value=self.sa_scores.rascore, unit="1")

        all_licenses = set(self.licenses)
        for off in self.offers:
            if off.license:
                all_licenses.add(off.license)
        for l_offs in self.leaf_offers.values():
            for lo in l_offs:
                if lo.license:
                    all_licenses.add(lo.license)
        governing_license = select_record_license(all_licenses, default="MIT")

        tags = [status_tag]
        if self.offers and not purchasable:
            tags.append("stock:all_offers_out_of_stock")
        if self.data_source != "unknown":
            tags.append(f"data_source:{self.data_source}")
        if self.inchikey_match != "none":
            tags.append(f"inchikey_match:{self.inchikey_match}")
        for lic in sorted(all_licenses):
            tags.append(f"license:{lic}")

        avail_hash = self.record_hash or self.content_hash()
        derived = list(dict.fromkeys([*self.provenance.derived_from, avail_hash]))
        prov = self.provenance.model_copy(update={"derived_from": derived})

        return AssayfoldRecord(
            record_id=make_assayfold_record_id(source_id, dataset_id, nid),
            source_id=source_id,
            dataset_id=dataset_id,
            native_id=nid,
            kind="molecule",
            payload=self.molecule,
            properties=props,
            blobs=[],
            license=governing_license,
            origin=self.origin,
            provenance=prov,
            tags=tags,
            schema_version=self.schema_version,
        ).sealed()


class RankResult(_Model):
    candidate: str
    name: Optional[str] = None
    canonical_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    status: Literal["buyable", "makeable", "unknown", "invalid"]
    n_steps: Optional[int] = None
    best_price: Optional[Quantity] = Field(
        default=None,
        description=(
            "Cheapest listed catalog price (normalized to '<currency>/g' when mass pack_size "
            "is known; retained as evidence even when all_offers_out_of_stock=True)."
        ),
    )
    est_cost: Optional[Quantity] = Field(
        default=None,
        description=(
            "Estimated stoichiometric leaf-only route cost per 1 g of target in '<currency>/g' "
            "(e.g. 'USD/g'); None for buyable, unknown, or invalid candidates."
        ),
    )
    est_cost_yield_adjusted: Optional[Quantity] = Field(
        default=None,
        description=(
            "Step-depth yield-compounded leaf-only route cost per 1 g of target in '<currency>/g'."
        ),
    )
    min_outlay_cost: Optional[Quantity] = Field(
        default=None,
        description=(
            "Minimum discrete pack MOQ cash outlay for the full target_scale_g order in "
            "absolute '<currency>' (e.g. 'USD')."
        ),
    )
    target_scale_g: float = Field(
        default=1.0,
        description="Target procurement or synthesis scale in grams used for min_outlay_cost.",
    )
    step_yield_assumed: float = 1.0
    cost_complete: bool = False
    all_offers_out_of_stock: bool = Field(
        default=False,
        description=(
            "True when >=1 catalog offer matched the molecule but every offer is explicitly "
            "flagged in_stock=False."
        ),
    )
    lead_time_days: Optional[int] = None
    moq_pack_g: Optional[float] = None
    precursor_smiles: list[str] = Field(default_factory=list)
    inchikey_match: InChIKeyMatchLiteral = "none"
    matched_offer_inchikeys: list[str] = Field(default_factory=list)
    stereoisomer_inchikeys_in_stock: list[str] = Field(default_factory=list)
    data_source: DataSourceLiteral = "unknown"
    sa_scores: Optional[SAScores] = None
    reasons: list[str] = Field(default_factory=list)
    evidence: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def _validate_rank_result(self) -> "RankResult":
        if (
            isinstance(self.target_scale_g, bool)
            or not math.isfinite(self.target_scale_g)
            or self.target_scale_g <= 0.0
        ):
            raise ValueError("RankResult.target_scale_g must be a finite number > 0")
        return self


PriceStatusLiteral = Literal[
    "verified",
    "leaf_lower_bound",
    "missing",
    "incomplete",
    "incompatible_unit",
    "not_applicable",
]


class ProposalDecision(_Model):
    candidate: str
    keep: bool
    over_budget: bool
    price_status: PriceStatusLiteral = "not_applicable"
    status: Literal["buyable", "makeable", "unknown", "invalid"]
    reasons: list[str] = Field(default_factory=list)
    result: RankResult


class FilterSummary(_Model):
    budget: Optional[float] = None
    budget_currency: str = "USD"
    require_complete_cost: bool = False
    max_outlay_budget: Optional[float] = None
    use_yield_adjusted_cost: bool = False
    kept: list[RankResult] = Field(default_factory=list)
    over_budget: list[RankResult] = Field(default_factory=list)
    rejected: list[RankResult] = Field(default_factory=list)
    decisions: list[ProposalDecision] = Field(default_factory=list)


def export_json_schema(out_dir: Path | str = "schema") -> list[Path]:
    """Write JSON Schema files for Source, Offer, Route, AvailabilityRecord, and RankResult."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    models = (
        ("source", Source),
        ("offer", Offer),
        ("route", Route),
        ("availability_record", AvailabilityRecord),
        ("rank_result", RankResult),
    )
    paths: list[Path] = []
    for name, model in models:
        p = out / f"{name}.schema.json"
        p.write_text(
            json.dumps(model.model_json_schema(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        paths.append(p)
    return paths
