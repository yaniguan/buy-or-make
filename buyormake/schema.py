"""Unified schema for buy-or-make (Pydantic v2, extra='forbid').

Field names follow the assayfold hub schema (Source, Quantity, Uncertainty,
AgentRef, SoftwareRef, Provenance, Molecule, origin, schema_version, record_hash)
so availability records can be exported to an open-data hub without renaming.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal, Optional

from pydantic import BaseModel, ConfigDict, Field

from .canonical import canonical_json, sha256_hex

SCHEMA_VERSION = "0.1.0"


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


class Quantity(_Model):
    value: Optional[float] = None
    unit: str = "1"
    raw_unit_str: Optional[str] = None
    uncertainty: Optional[Uncertainty] = None


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
    license_url: Optional[str] = None
    license_evidence: str
    license_quote: Optional[str] = None
    license_checked: str
    commercial_use_ok: bool
    redistribution_ok: bool
    attribution: str
    modalities: list[str] = Field(default_factory=list)
    kinds: list[str] = Field(default_factory=list)
    access: Literal["api", "download", "huggingface", "github", "builtin"]
    size_full: Optional[str] = None
    status: Literal["integrated", "fetch_only", "not_used", "planned"] = "integrated"
    notes: Optional[str] = None


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


class RouteStep(_Model):
    product: str
    reactants: list[str]
    reaction_smiles: Optional[str] = None
    template_or_reference: str
    source: str


class Route(_Model):
    target_inchikey: str
    n_steps: int
    steps: list[RouteStep]
    leaves: list[str]
    leaves_purchasable: list[bool]
    method: str
    score: Optional[float] = None
    est_cost: Optional[Quantity] = None
    est_cost_method: str
    est_cost_complete: bool


class SAScores(_Model):
    sascore: Optional[float] = None
    scscore: Optional[float] = None
    rascore: Optional[float] = None


OriginLiteral = Literal["real", "simulated", "model_generated", "synthetic_demo"]


class AvailabilityRecord(_Model):
    molecule: Molecule
    offers: list[Offer] = Field(default_factory=list)
    routes: list[Route] = Field(default_factory=list)
    sa_scores: SAScores = Field(default_factory=SAScores)
    origin: OriginLiteral = "real"
    provenance: Provenance = Field(default_factory=Provenance)
    licenses: list[str] = Field(default_factory=list)
    schema_version: str = SCHEMA_VERSION
    record_hash: Optional[str] = None

    def content_hash(self) -> str:
        """SHA-256 over canonical JSON. Excludes record_hash and provenance.recorded_utc."""
        d = self.model_dump(mode="json", exclude={"record_hash"})
        d["provenance"].pop("recorded_utc", None)
        return sha256_hex(canonical_json(d))

    def sealed(self) -> "AvailabilityRecord":
        """Return a copy with record_hash filled in."""
        return self.model_copy(update={"record_hash": self.content_hash()})


class RankResult(_Model):
    candidate: str
    canonical_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    status: Literal["buyable", "makeable", "unknown", "invalid"]
    n_steps: Optional[int] = None
    best_price: Optional[Quantity] = None
    est_cost: Optional[Quantity] = None
    cost_complete: bool = False
    sa_scores: Optional[SAScores] = None
    reasons: list[str] = Field(default_factory=list)
    evidence: list[str] = Field(default_factory=list)


class ProposalDecision(_Model):
    candidate: str
    keep: bool
    over_budget: bool
    status: Literal["buyable", "makeable", "unknown", "invalid"]
    reasons: list[str] = Field(default_factory=list)
    result: RankResult


class FilterSummary(_Model):
    budget: Optional[float] = None
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
        p.write_text(json.dumps(model.model_json_schema(), indent=2, sort_keys=True) + "\n")
        paths.append(p)
    return paths
