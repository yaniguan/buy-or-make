"""Frontier multi-step retrosynthesis: Chen et al. (ICML 2020) `Retro*` best-first AND-OR search
with `SCScore`/`SAscore` value heuristics, extended atom-conserved medicinal chemistry SMARTS,
and exact discrete pack-size MOQ outlay minimization across shared multi-branch building blocks.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping, Optional

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .canonical import NormalizedMol, normalize_smiles
from .cost import (
    _molar_mass_ratio,
    compute_leaf_depths,
    estimate_route_cost_bounds,
    offer_pack_grams,
    route_cost_method,
)
from .routes import (
    RETRO_TEMPLATE_DELTAS,
    RETRO_TEMPLATES,
    UNIMOLECULAR_RING_OPENING_TEMPLATES,
    ReactionLookupBackend,
    TemplateFallbackBackend,
    _has_chemoselectivity_conflict,
    _resolve_leaf_offers,
)
from .schema import Offer, Quantity, Route, RouteStep
from .scores import DEFAULT_SCSCORE_WEIGHTS, compute_sascore, compute_scscore

OfferLookupFn = Callable[[str, str], list[Offer]]

# Extended atom-conserved medicinal chemistry SMARTS templates complementing RETRO_TEMPLATES.
# Every bimolecular template enforces `!@` on the cleaved bond; unimolecular FGIs / hydrolytic
# ring-openings are registered explicitly with exact heavy-atom conservation deltas.
EXTENDED_RETRO_STAR_TEMPLATES: list[tuple[str, str]] = [
    *RETRO_TEMPLATES,
    (
        "suzuki_biaryl_cl",
        "[c:1]-!@[c:2]>>[c:1]Cl.[c:2]B(O)O",
    ),
    (
        "suzuki_biaryl_i",
        "[c:1]-!@[c:2]>>[c:1]I.[c:2]B(O)O",
    ),
    (
        "buchwald_hartwig_cn_cl",
        "[c:1]!@[NH1:2][c:3]>>[c:1]Cl.[NH2:2][c:3]",
    ),
    (
        "snar_heteroaryl_cl_tertiary",
        "[c;$(c1ncncn1),$(c1ncccn1),$(c1ccncn1),$(c1ncncc1),$(c1ccccn1):1]!@[N:2]([#6:3])[#6:4]>>[c:1]Cl.[NH1:2]([#6:3])[#6:4]",
    ),
    (
        "snar_heteroaryl_cl_secondary",
        "[c;$(c1ncncn1),$(c1ncccn1),$(c1ccncn1),$(c1ncncc1),$(c1ccccn1):1]!@[NH1:2][#6:3]>>[c:1]Cl.[NH2:2][#6:3]",
    ),
    (
        "acyl_chloride_amide_secondary",
        "[#6:1][C:2](=[O:3])!@[NH1:4][#6:5]>>[#6:1][C:2](=[O:3])Cl.[NH2:4][#6:5]",
    ),
    (
        "acyl_chloride_amide_tertiary",
        "[#6:1][C:2](=[O:3])!@[N:4]([#6:5])[#6:6]>>[#6:1][C:2](=[O:3])Cl.[NH1:4]([#6:5])[#6:6]",
    ),
    (
        "n_alkylation_secondary_br",
        "[#6:1][NH1;!$(NC=[O,S]);!$(NS(=O)=O):2]!@[CX4;H3,H2;!$(C([#6])([#6])[#6]):3]>>[#6:1][NH2:2].[CX4:3]Br",
    ),
    (
        "n_alkylation_tertiary_br",
        "[N+0;H0;!a;!$(NC=[O,S]);!$(NS(=O)=O);!$([N+](=O)[O-]):1]([#6:4])([#6:5])!@[CX4;H3,H2;!$(C([#6])([#6])[#6]):2]>>[NH1:1]([#6:4])[#6:5].[CX4:2]Br",
    ),
    (
        "n_alkylation_tertiary_cl",
        "[N+0;H0;!a;!$(NC=[O,S]);!$(NS(=O)=O);!$([N+](=O)[O-]):1]([#6:4])([#6:5])!@[CX4;H3,H2;!$(C([#6])([#6])[#6]):2]>>[NH1:1]([#6:4])[#6:5].[CX4:2]Cl",
    ),
    (
        "n_heterocycle_alkylation_br",
        "[nH0;+0;r5,r6:1]!@[CX4;H3,H2;!$(C([#6])([#6])[#6]):2]>>[nH1:1].[CX4:2]Br",
    ),
    (
        "williamson_aryl_ether_cl",
        "[c:1][O:2]!@[CX4;H3,H2;!$(C([#6])([#6])[#6]):3]>>[c:1][OH1:2].[CX4:3]Cl",
    ),
    (
        "thioether_alkylation_br",
        "[#6:1][SX2:2]!@[CX4;H3,H2;!$(C([#6])([#6])[#6]):3]>>[#6:1][SH1:2].[CX4:3]Br",
    ),
    (
        "sonogashira_alkynyl",
        "[c:1]-!@[C:2]#[C:3]>>[c:1]Br.[CH1:2]#[C:3]",
    ),
    (
        "n_boc_protection_secondary",
        "[#6:1][NH1:2]!@[C:3](=[O:4])[O:5][C:6]([CH3:7])([CH3:8])[CH3:9]>>[#6:1][NH2:2].CC(C)(C)OC(=O)[O:5][C:3](=[O:4])O[C:6]([CH3:7])([CH3:8])[CH3:9]",
    ),
    (
        "n_boc_protection_tertiary",
        "[N+0;H0;!a:1]([#6:10])([#6:11])!@[C:3](=[O:4])[O:5][C:6]([CH3:7])([CH3:8])[CH3:9]>>[NH1:1]([#6:10])[#6:11].CC(C)(C)OC(=O)[O:5][C:3](=[O:4])O[C:6]([CH3:7])([CH3:8])[CH3:9]",
    ),
    (
        "n_boc_protection_heteroaryl",
        "[n+0;H0;r5,r6:1]!@[C:3](=[O:4])[O:5][C:6]([CH3:7])([CH3:8])[CH3:9]>>[nH1:1].CC(C)(C)OC(=O)[O:5][C:3](=[O:4])O[C:6]([CH3:7])([CH3:8])[CH3:9]",
    ),
    (
        "nitro_reduction_aryl",
        "[c:1][NH2:2]>>[c:1][N+:2](=O)[O-]",
    ),
    (
        "ester_hydrolysis_methyl",
        "[#6:1][C:2](=[O:3])[OH1:4]>>[#6:1][C:2](=[O:3])[O:4]C",
    ),
    (
        "ester_hydrolysis_ethyl",
        "[#6:1][C:2](=[O:3])[OH1:4]>>[#6:1][C:2](=[O:3])[O:4]CC",
    ),
    (
        "ester_hydrolysis_tbu",
        "[#6:1][C:2](=[O:3])[OH1:4]>>[#6:1][C:2](=[O:3])[O:4]C(C)(C)C",
    ),
    (
        "ketone_reduction_secondary_alcohol",
        "[#6:1][CH1:2]([OH1:3])[#6:4]>>[#6:1][C:2](=[O:3])[#6:4]",
    ),
    (
        "aldehyde_reduction_primary_alcohol",
        "[#6:1][CH2:2][OH1:3]>>[#6:1][CH1:2]=[O:3]",
    ),
    (
        "boc_deprotection_primary_amine",
        "[#6:1][NH2:2]>>[#6:1][NH1:2]C(=O)OC(C)(C)C",
    ),
    (
        "boc_deprotection_secondary_amine",
        "[#6:1][NH1;!$(NC=O);!$(NS=O):2][#6:3]>>[#6:1][N:2]([#6:3])C(=O)OC(C)(C)C",
    ),
    (
        "cbz_deprotection_primary_amine",
        "[#6:1][NH2:2]>>[#6:1][NH1:2]C(=O)OCc1ccccc1",
    ),
    (
        "cbz_deprotection_secondary_amine",
        "[#6:1][NH1;!$(NC=O);!$(NS=O):2][#6:3]>>[#6:1][N:2]([#6:3])C(=O)OCc1ccccc1",
    ),
    (
        "bn_ether_deprotection_phenol",
        "[c:1][OH1:2]>>[c:1][O:2]Cc1ccccc1",
    ),
    (
        "bn_ether_deprotection_alcohol",
        "[CX4:1][OH1:2]>>[CX4:1][O:2]Cc1ccccc1",
    ),
    (
        "ester_hydrolysis_benzyl",
        "[#6:1][C:2](=[O:3])[OH1:4]>>[#6:1][C:2](=[O:3])[O:4]Cc1ccccc1",
    ),
    (
        "nitrile_reduction_primary_amine",
        "[#6:1][CH2:2][NH2:3]>>[#6:1][C:2]#[N:3]",
    ),
    (
        "ester_reduction_primary_alcohol_me",
        "[#6:1][CH2:2][OH1:3]>>[#6:1][C:2](=[O:3])OC",
    ),
    (
        "snar_activated_aryl_f_secondary",
        "[c;$(c1ccc([N+](=O)[O-])cc1),$(c1ccc(C#N)cc1),$(c1ccccn1),$(c1ncccn1):1]!@[NH1:2][#6:3]>>[c:1]F.[NH2:2][#6:3]",
    ),
    (
        "snar_activated_aryl_f_tertiary",
        "[c;$(c1ccc([N+](=O)[O-])cc1),$(c1ccc(C#N)cc1),$(c1ccccn1),$(c1ncccn1):1]!@[N:2]([#6:3])[#6:4]>>[c:1]F.[NH1:2]([#6:3])[#6:4]",
    ),
]

RETRO_STAR_TEMPLATE_DELTAS: dict[str, dict[str, int]] = {
    **RETRO_TEMPLATE_DELTAS,
    "suzuki_biaryl_cl": {"Cl": 1, "B": 1, "O": 2},
    "suzuki_biaryl_i": {"I": 1, "B": 1, "O": 2},
    "buchwald_hartwig_cn_cl": {"Cl": 1},
    "snar_heteroaryl_cl_tertiary": {"Cl": 1},
    "snar_heteroaryl_cl_secondary": {"Cl": 1},
    "acyl_chloride_amide_secondary": {"Cl": 1},
    "acyl_chloride_amide_tertiary": {"Cl": 1},
    "n_alkylation_secondary_br": {"Br": 1},
    "n_alkylation_tertiary_br": {"Br": 1},
    "n_alkylation_tertiary_cl": {"Cl": 1},
    "n_heterocycle_alkylation_br": {"Br": 1},
    "williamson_aryl_ether_cl": {"Cl": 1},
    "thioether_alkylation_br": {"Br": 1},
    "sonogashira_alkynyl": {"Br": 1},
    "n_boc_protection_secondary": {"C": 5, "O": 3},
    "n_boc_protection_tertiary": {"C": 5, "O": 3},
    "n_boc_protection_heteroaryl": {"C": 5, "O": 3},
    "nitro_reduction_aryl": {"O": 2},
    "ester_hydrolysis_methyl": {"C": 1},
    "ester_hydrolysis_ethyl": {"C": 2},
    "ester_hydrolysis_tbu": {"C": 4},
    "ketone_reduction_secondary_alcohol": {},
    "aldehyde_reduction_primary_alcohol": {},
    "boc_deprotection_primary_amine": {"C": 5, "O": 2},
    "boc_deprotection_secondary_amine": {"C": 5, "O": 2},
    "cbz_deprotection_primary_amine": {"C": 8, "O": 2},
    "cbz_deprotection_secondary_amine": {"C": 8, "O": 2},
    "bn_ether_deprotection_phenol": {"C": 7},
    "bn_ether_deprotection_alcohol": {"C": 7},
    "ester_hydrolysis_benzyl": {"C": 7},
    "nitrile_reduction_primary_amine": {},
    "ester_reduction_primary_alcohol_me": {"C": 1, "O": 1},
    "snar_activated_aryl_f_secondary": {"F": 1},
    "snar_activated_aryl_f_tertiary": {"F": 1},
}

UNIMOLECULAR_RETRO_STAR_TEMPLATES: frozenset[str] = frozenset(
    {
        *UNIMOLECULAR_RING_OPENING_TEMPLATES,
        "nitro_reduction_aryl",
        "ester_hydrolysis_methyl",
        "ester_hydrolysis_ethyl",
        "ester_hydrolysis_tbu",
        "ester_hydrolysis_benzyl",
        "ketone_reduction_secondary_alcohol",
        "aldehyde_reduction_primary_alcohol",
        "nitrile_reduction_primary_amine",
        "ester_reduction_primary_alcohol_me",
        "boc_deprotection_primary_amine",
        "boc_deprotection_secondary_amine",
        "cbz_deprotection_primary_amine",
        "cbz_deprotection_secondary_amine",
        "bn_ether_deprotection_phenol",
        "bn_ether_deprotection_alcohol",
    }
)

_SIZE_INCREASING_FGI_TEMPLATES: frozenset[str] = frozenset(
    {
        "nitro_reduction_aryl",
        "ester_hydrolysis_methyl",
        "ester_hydrolysis_ethyl",
        "ester_hydrolysis_tbu",
        "ester_hydrolysis_benzyl",
        "ester_reduction_primary_alcohol_me",
        "boc_deprotection_primary_amine",
        "boc_deprotection_secondary_amine",
        "cbz_deprotection_primary_amine",
        "cbz_deprotection_secondary_amine",
        "bn_ether_deprotection_phenol",
        "bn_ether_deprotection_alcohol",
    }
)

_PAT_UNPROTECTED_AMINE = Chem.MolFromSmarts(
    "[N;!a;H2,H1;!$(NC(=[O,S]));!$(NS(=O)=O);!$(N=[C,N,O]);!$(N[N,O])]"
)
_PAT_UNPROTECTED_OH = Chem.MolFromSmarts(
    "[OH1;!$(OC(=[O,S]));!$(OS(=O));!$(OB);!$(OP);!$(ON)]"
)
_PAT_ARYL_HAL_ANY = Chem.MolFromSmarts("[c,CX3][Cl,Br,I]")
_PAT_BORONIC_ACID = Chem.MolFromSmarts("[#6]B([OH])[OH]")


@dataclass(frozen=True)
class _CompiledExtTemplate:
    name: str
    smarts: str
    rxn: rdChemReactions.ChemicalReaction
    heavy_delta: dict[str, int]
    base_cost: float


_COMPILED_EXT_TEMPLATES: Optional[list[_CompiledExtTemplate]] = None

_TEMPLATE_BASE_COSTS: dict[str, float] = {
    "amide_coupling_secondary": 1.10,
    "amide_coupling_tertiary": 1.12,
    "acyl_chloride_amide_secondary": 1.18,
    "acyl_chloride_amide_tertiary": 1.20,
    "esterification": 1.15,
    "sulfonamide_secondary": 1.15,
    "sulfonamide_tertiary": 1.18,
    "urea_formation": 1.20,
    "carbamate_formation": 1.22,
    "suzuki_biaryl": 1.12,
    "suzuki_biaryl_cl": 1.22,
    "suzuki_biaryl_i": 1.24,
    "suzuki_aryl_alkenyl": 1.22,
    "buchwald_hartwig_cn": 1.24,
    "buchwald_hartwig_cn_tertiary": 1.26,
    "buchwald_hartwig_cn_cl": 1.30,
    "snar_heteroaryl_cl_tertiary": 1.16,
    "snar_heteroaryl_cl_secondary": 1.16,
    "williamson_aryl_ether": 1.18,
    "williamson_aryl_ether_cl": 1.24,
    "williamson_alkyl_ether": 1.26,
    "n_alkylation_secondary_br": 1.24,
    "n_alkylation_tertiary_br": 1.20,
    "n_alkylation_tertiary_cl": 1.26,
    "n_heterocycle_alkylation_br": 1.22,
    "thioether_alkylation_br": 1.24,
    "reductive_amination_secondary": 1.20,
    "reductive_amination_tertiary": 1.22,
    "sonogashira_alkynyl": 1.28,
    "n_boc_protection_secondary": 1.14,
    "n_boc_protection_tertiary": 1.14,
    "n_boc_protection_heteroaryl": 1.16,
    "lactone_ring_opening": 1.25,
    "lactam_ring_opening": 1.28,
    "nitro_reduction_aryl": 1.20,
    "ester_hydrolysis_methyl": 1.22,
    "ester_hydrolysis_ethyl": 1.22,
    "ester_hydrolysis_tbu": 1.24,
    "ketone_reduction_secondary_alcohol": 1.26,
    "aldehyde_reduction_primary_alcohol": 1.26,
    "boc_deprotection_primary_amine": 1.20,
    "boc_deprotection_secondary_amine": 1.20,
    "cbz_deprotection_primary_amine": 1.24,
    "bn_ether_deprotection_phenol": 1.26,
}


def _get_compiled_ext_templates() -> list[_CompiledExtTemplate]:
    global _COMPILED_EXT_TEMPLATES
    if _COMPILED_EXT_TEMPLATES is None:
        compiled: list[_CompiledExtTemplate] = []
        for name, smarts in EXTENDED_RETRO_STAR_TEMPLATES:
            rxn = rdChemReactions.ReactionFromSmarts(smarts)
            compiled.append(
                _CompiledExtTemplate(
                    name=name,
                    smarts=smarts,
                    rxn=rxn,
                    heavy_delta=dict(RETRO_STAR_TEMPLATE_DELTAS[name]),
                    base_cost=_TEMPLATE_BASE_COSTS.get(name, 1.30),
                )
            )
        _COMPILED_EXT_TEMPLATES = compiled
    return _COMPILED_EXT_TEMPLATES


def _heavy_element_counts(mol: Chem.Mol) -> Counter[str]:
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def _has_extended_chemoselectivity_conflict(tpl_name: str, frag_mols: list[Chem.Mol]) -> bool:
    """Enforce functional-group chemoselectivity and AB-homopolymerization guards."""
    if _has_chemoselectivity_conflict(tpl_name, frag_mols):
        return True
    if len(frag_mols) < 2:
        return False
    f0, f1 = frag_mols[0], frag_mols[1]

    if tpl_name in {"acyl_chloride_amide_secondary", "acyl_chloride_amide_tertiary"}:
        if f0.HasSubstructMatch(_PAT_UNPROTECTED_AMINE) or f0.HasSubstructMatch(_PAT_UNPROTECTED_OH):
            return True
        if len(f1.GetSubstructMatches(_PAT_UNPROTECTED_AMINE)) >= 2:
            return True
        return False

    if tpl_name in {"suzuki_biaryl_cl", "suzuki_biaryl_i"}:
        if f0.HasSubstructMatch(_PAT_BORONIC_ACID):
            return True
        if f1.HasSubstructMatch(_PAT_ARYL_HAL_ANY):
            return True
        if len(f1.GetSubstructMatches(_PAT_BORONIC_ACID)) >= 2:
            return True
        if len(f0.GetSubstructMatches(_PAT_ARYL_HAL_ANY)) >= 2:
            return True
        return False

    if tpl_name in {
        "buchwald_hartwig_cn_cl",
        "snar_heteroaryl_cl_tertiary",
        "snar_heteroaryl_cl_secondary",
    }:
        if f0.HasSubstructMatch(_PAT_UNPROTECTED_AMINE):
            return True
        return False

    if tpl_name in {
        "n_alkylation_secondary_br",
        "n_alkylation_tertiary_br",
        "n_alkylation_tertiary_cl",
        "n_heterocycle_alkylation_br",
        "williamson_aryl_ether_cl",
        "thioether_alkylation_br",
    }:
        if f1.HasSubstructMatch(_PAT_UNPROTECTED_AMINE) or f1.HasSubstructMatch(_PAT_UNPROTECTED_OH):
            return True
        ring_info = f1.GetRingInfo()
        for atom in f1.GetAtoms():
            if atom.GetAtomicNum() in {17, 35}:
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6:
                        if nbr.GetTotalNumHs() == 0:
                            return True
                        if ring_info.NumAtomRings(nbr.GetIdx()) >= 2:
                            return True
        return False

    return False


# ---------- Exact Discrete Pack-Size MOQ Outlay Optimizer ----------


@dataclass(frozen=True)
class PackAllocation:
    """Optimal discrete SKU pack purchase entry for a single leaf building block."""

    leaf_smiles: str
    leaf_inchikey: str
    required_grams: float
    purchased_grams: float
    outlay_usd: float
    vendor: str
    catalog_id: Optional[str]
    pack_grams: float
    pack_count: int
    unit_pack_price_usd: float


@dataclass(frozen=True)
class MOQOptimizationResult:
    """Exact shared-synthon discrete pack-size MOQ outlay across a multi-step route."""

    pooled_outlay: Optional[Quantity]
    unpooled_outlay: Optional[Quantity]
    pooling_savings_usd: float
    complete: bool
    shared_leaf_count: int
    allocations: list[PackAllocation] = field(default_factory=list)


def optimize_discrete_pack_order(
    offers: list[Offer],
    *,
    required_grams: float,
    currency: str = "USD",
    target_mol_weight: Optional[float] = None,
    target_inchikey: Optional[str] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
) -> tuple[Optional[float], list[tuple[Offer, float, int]]]:
    """Solve exact integer pack-cover minimizing cash outlay `sum(n_i * price_i)` for `required_grams`.

    Unlike naive single-SKU ceiling division, evaluates both single-SKU and mixed-SKU pack
    combinations (e.g., `1 x 1 g + 1 x 250 mg` vs `2 x 1 g`) via LP-pruned branch-and-bound.
    Returns `(min_outlay_in_currency, [(offer, effective_pack_g, count), ...])`.
    """
    if required_grams <= 0 or not math.isfinite(required_grams):
        raise ValueError("required_grams must be a finite number > 0")

    # Convert each priced mass SKU into effective target-equivalent grams `eff_pack_g`
    raw_skus: list[tuple[float, float, Offer]] = []  # (eff_pack_g, price, offer)
    for off in offers:
        if off.price is None or off.currency != currency:
            continue
        pack_g = offer_pack_grams(off)
        if pack_g is None or pack_g <= 0:
            continue
        mw_ratio, _ = _molar_mass_ratio(
            off.inchikey, target_inchikey, target_mol_weight, offer_mol_weights
        )
        eff_pack_g = pack_g / max(mw_ratio, 1e-12)
        if eff_pack_g <= 0:
            continue
        raw_skus.append((eff_pack_g, float(off.price), off))

    if not raw_skus:
        return None, []

    # Deduplicate identical effective pack sizes by keeping cheapest price
    by_size_ug: dict[int, tuple[float, float, Offer]] = {}
    for eff_g, price, off in raw_skus:
        key_ug = max(1, int(round(eff_g * 1_000_000.0)))
        prev = by_size_ug.get(key_ug)
        if prev is None or price < prev[1]:
            by_size_ug[key_ug] = (eff_g, price, off)

    # Sort SKUs by unit price (price / eff_pack_g) ascending for tight LP lower-bound pruning
    skus = sorted(by_size_ug.values(), key=lambda item: (item[1] / item[0], -item[0]))
    skus = skus[:16]  # Bound distinct non-identical pack tiers per chemical

    # Seed incumbent with best single-SKU pack cover
    best_cost = math.inf
    best_plan: list[tuple[Offer, float, int]] = []
    for eff_g, price, off in skus:
        n_single = max(1, int(math.ceil(required_grams / eff_g - 1e-12)))
        c_single = n_single * price
        if c_single < best_cost - 1e-12:
            best_cost = c_single
            best_plan = [(off, eff_g, n_single)]

    if len(skus) == 1 or best_cost == 0.0:
        return round(best_cost, 4), best_plan

    # Suffix minimum unit-price per gram for admissible LP relaxation pruning
    n_skus = len(skus)
    suffix_min_ppg = [0.0] * n_skus
    running = math.inf
    for i in range(n_skus - 1, -1, -1):
        running = min(running, skus[i][1] / skus[i][0])
        suffix_min_ppg[i] = running

    cur_counts = [0] * n_skus

    def _dfs(idx: int, rem_g: float, acc_cost: float) -> None:
        nonlocal best_cost, best_plan
        if rem_g <= 1e-12:
            if acc_cost < best_cost - 1e-12:
                best_cost = acc_cost
                best_plan = [
                    (skus[j][2], skus[j][0], cur_counts[j])
                    for j in range(n_skus)
                    if cur_counts[j] > 0
                ]
            return
        if idx >= n_skus:
            return
        # Admissible fractional lower bound
        if acc_cost + rem_g * suffix_min_ppg[idx] >= best_cost - 1e-12:
            return

        eff_g, price, _ = skus[idx]
        if idx == n_skus - 1:
            n_last = max(1, int(math.ceil(rem_g / eff_g - 1e-12)))
            cand_cost = acc_cost + n_last * price
            if cand_cost < best_cost - 1e-12:
                cur_counts[idx] = n_last
                best_cost = cand_cost
                best_plan = [
                    (skus[j][2], skus[j][0], cur_counts[j])
                    for j in range(n_skus)
                    if cur_counts[j] > 0
                ]
                cur_counts[idx] = 0
            return

        max_n = int(math.ceil(rem_g / eff_g - 1e-12))
        if price > 0:
            max_n = min(max_n, int((best_cost - acc_cost) // price) + 1)
        # Cap branching on tiny milligram vials when required_grams is huge
        if max_n > 64:
            candidate_ns = [max_n, 0]
        else:
            candidate_ns = list(range(max_n, -1, -1))

        for n_take in candidate_ns:
            next_cost = acc_cost + n_take * price
            if next_cost >= best_cost - 1e-12:
                continue
            cur_counts[idx] = n_take
            _dfs(idx + 1, rem_g - n_take * eff_g, next_cost)
        cur_counts[idx] = 0

    _dfs(0, required_grams, 0.0)
    return round(best_cost, 4), best_plan


def optimize_route_moq_outlay(
    target_smiles: str,
    leaves: list[str],
    leaf_offers: Mapping[str, list[Offer]],
    *,
    steps: Optional[list[RouteStep]] = None,
    currency: str = "USD",
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
    target_scale_g: float = 1.0,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
) -> MOQOptimizationResult:
    """Minimize discrete pack-size MOQ cash outlay across pooled multi-branch building blocks."""
    target_nm = normalize_smiles(target_smiles)
    if not leaves or target_nm is None or target_nm.mol_weight <= 0:
        return MOQOptimizationResult(
            pooled_outlay=None,
            unpooled_outlay=None,
            pooling_savings_usd=0.0,
            complete=False,
            shared_leaf_count=0,
        )

    depth_pairs = compute_leaf_depths(target_smiles, leaves, steps)
    depths_by_can: dict[str, list[int]] = {}
    nm_by_can: dict[str, NormalizedMol] = {}
    for raw_leaf, depth in depth_pairs:
        lnm = normalize_smiles(raw_leaf)
        if lnm is None:
            return MOQOptimizationResult(
                pooled_outlay=None,
                unpooled_outlay=None,
                pooling_savings_usd=0.0,
                complete=False,
                shared_leaf_count=0,
            )
        depths_by_can.setdefault(lnm.canonical_smiles, []).append(depth)
        nm_by_can[lnm.canonical_smiles] = lnm

    total_pooled = 0.0
    total_unpooled = 0.0
    priced_unique = 0
    shared_leaves = 0
    allocations: list[PackAllocation] = []

    for can_smi, depths in depths_by_can.items():
        if len(depths) >= 2:
            shared_leaves += 1
        lnm = nm_by_can[can_smi]
        offs = list(leaf_offers.get(lnm.inchikey, []))
        mw_frac = lnm.mol_weight / target_nm.mol_weight
        branch_masses = [
            max(1e-9, target_scale_g * (stoich_excess / (default_step_yield**d)) * mw_frac)
            for d in depths
        ]
        pooled_req_g = sum(branch_masses)

        pooled_cost, plan = optimize_discrete_pack_order(
            offs,
            required_grams=pooled_req_g,
            currency=currency,
            target_mol_weight=lnm.mol_weight,
            target_inchikey=lnm.inchikey,
            offer_mol_weights=offer_mol_weights,
        )
        if pooled_cost is None:
            continue

        unpooled_leaf_cost = 0.0
        for b_mass in branch_masses:
            b_cost, _ = optimize_discrete_pack_order(
                offs,
                required_grams=b_mass,
                currency=currency,
                target_mol_weight=lnm.mol_weight,
                target_inchikey=lnm.inchikey,
                offer_mol_weights=offer_mol_weights,
            )
            if b_cost is not None:
                unpooled_leaf_cost += b_cost

        total_pooled += pooled_cost
        total_unpooled += max(pooled_cost, unpooled_leaf_cost)
        priced_unique += 1
        for sel_off, eff_pack_g, count in plan:
            allocations.append(
                PackAllocation(
                    leaf_smiles=can_smi,
                    leaf_inchikey=lnm.inchikey,
                    required_grams=round(pooled_req_g, 6),
                    purchased_grams=round(eff_pack_g * count, 6),
                    outlay_usd=round(float(sel_off.price or 0.0) * count, 4),
                    vendor=sel_off.vendor,
                    catalog_id=sel_off.catalog_id,
                    pack_grams=round(eff_pack_g, 6),
                    pack_count=count,
                    unit_pack_price_usd=round(float(sel_off.price or 0.0), 4),
                )
            )

    if priced_unique == 0:
        return MOQOptimizationResult(
            pooled_outlay=None,
            unpooled_outlay=None,
            pooling_savings_usd=0.0,
            complete=False,
            shared_leaf_count=shared_leaves,
        )

    complete = priced_unique == len(depths_by_can)
    savings = max(0.0, round(total_unpooled - total_pooled, 4))
    note_pool = (
        f"; saved {savings:.4f} {currency} by pooling {shared_leaves} shared leaf(ves)"
        if savings > 0
        else ""
    )
    pooled_q = Quantity(
        value=round(total_pooled, 4),
        unit=currency,
        raw_unit_str=(
            f"{round(total_pooled, 4)} {currency} "
            f"(discrete pack MOQ outlay for {target_scale_g:g} g target{note_pool})"
        ),
    )
    unpooled_q = Quantity(
        value=round(total_unpooled, 4),
        unit=currency,
        raw_unit_str=f"{round(total_unpooled, 4)} {currency} (unpooled per-branch MOQ outlay)",
    )
    return MOQOptimizationResult(
        pooled_outlay=pooled_q,
        unpooled_outlay=unpooled_q,
        pooling_savings_usd=savings,
        complete=complete,
        shared_leaf_count=shared_leaves,
        allocations=allocations,
    )


# ---------- Chen et al. (ICML 2020) Retro* Best-First AND-OR Tree Search ----------


@dataclass
class _OrNode:
    node_id: int
    smiles: str
    inchikey: str
    mol: Optional[Chem.Mol]
    depth: int
    ancestors: frozenset[str]
    parent: Optional["_AndNode"]
    in_stock: bool
    v0: float
    leaf_cost_bonus: float = 0.0
    expanded: bool = False
    is_dead: bool = False
    succ: bool = False
    v_self: float = 0.0
    min_steps: int = 1
    children: list["_AndNode"] = field(default_factory=list)


@dataclass
class _AndNode:
    rxn_id: int
    step: RouteStep
    cost: float
    parent: _OrNode
    children: list[_OrNode] = field(default_factory=list)
    succ: bool = False
    v_self: float = math.inf
    min_steps: int = 1


class RetroStarPlanner:
    """Chen et al. (ICML 2020) `Retro*` best-first AND-OR retrosynthesis planner.

    Maintains target-rooted synthesis cost numbers `V_t(m | T)` and admissible minimum
    reaction-step counts across an explicit cycle-free AND-OR tree up to `max_steps` (default
    supporting 1..6+ steps, both linear and multi-branch convergent routes).
    """

    def __init__(
        self,
        *,
        max_expansions: int = 90,
        max_branch_per_node: int = 20,
        post_solve_expansions: int = 18,
        top_k_routes: int = 5,
        use_extended_templates: bool = True,
        scscore_weights: Path | str = DEFAULT_SCSCORE_WEIGHTS,
    ) -> None:
        self.max_expansions = max(1, int(max_expansions))
        self.max_branch_per_node = max(1, int(max_branch_per_node))
        self.post_solve_expansions = max(0, int(post_solve_expansions))
        self.top_k_routes = max(1, int(top_k_routes))
        self.use_extended_templates = bool(use_extended_templates)
        self.scscore_weights = Path(scscore_weights)
        self._ext_templates = _get_compiled_ext_templates() if self.use_extended_templates else []

    def expand_extended_smarts(
        self,
        canonical_smiles: str,
        *,
        mol: Optional[Chem.Mol] = None,
        is_stock_fn: Optional[Callable[[str], bool]] = None,
        allow_unstocked_fgi: bool = False,
    ) -> list[tuple[float, RouteStep]]:
        """Generate atom-conserved, chemoselectivity-filtered 1-step SMARTS disconnections."""
        parsed = mol if mol is not None else Chem.MolFromSmiles(canonical_smiles)
        if parsed is None:
            return []
        parent_heavy = parsed.GetNumHeavyAtoms()
        parent_counts = _heavy_element_counts(parsed)
        candidates: list[tuple[tuple[int, float, int], float, RouteStep]] = []
        seen_sets: set[tuple[str, ...]] = set()

        for tpl in self._ext_templates:
            try:
                outcomes = tpl.rxn.RunReactants((parsed,))
            except Exception:
                continue
            expected_counts = Counter(parent_counts)
            for sym, delta in tpl.heavy_delta.items():
                expected_counts[sym] += delta
                if expected_counts[sym] <= 0:
                    del expected_counts[sym]

            is_size_inc_fgi = tpl.name in _SIZE_INCREASING_FGI_TEMPLATES
            max_allowed_frag = parent_heavy + 9 if is_size_inc_fgi else parent_heavy + 5

            for prod_tuple in outcomes:
                reactants: list[str] = []
                frag_mols: list[Chem.Mol] = []
                frag_sizes: list[int] = []
                frag_counts: Counter[str] = Counter()
                ok = True
                for frag in prod_tuple:
                    try:
                        Chem.SanitizeMol(frag)
                        smi = Chem.MolToSmiles(frag, isomericSmiles=True, canonical=True)
                    except Exception:
                        ok = False
                        break
                    if "." in smi:
                        ok = False
                        break
                    nm = normalize_smiles(smi)
                    if nm is None or nm.canonical_smiles == canonical_smiles or "." in nm.canonical_smiles:
                        ok = False
                        break
                    if nm.num_heavy_atoms < 2 or nm.mol is None:
                        ok = False
                        break
                    frag_mols.append(nm.mol)
                    frag_sizes.append(nm.num_heavy_atoms)
                    frag_counts.update(_heavy_element_counts(nm.mol))
                    reactants.append(nm.canonical_smiles)

                if not ok or not reactants:
                    continue
                if frag_counts != expected_counts:
                    continue
                if any(sz > max_allowed_frag for sz in frag_sizes):
                    continue
                if _has_extended_chemoselectivity_conflict(tpl.name, frag_mols):
                    continue
                if is_size_inc_fgi and not allow_unstocked_fgi and is_stock_fn is not None:
                    if not all(is_stock_fn(r) for r in reactants):
                        continue

                key = tuple(sorted(reactants))
                if key in seen_sets:
                    continue
                seen_sets.add(key)

                stock_misses = (
                    sum(0 if is_stock_fn(r) else 1 for r in reactants)
                    if is_stock_fn is not None
                    else len(reactants)
                )
                step_cost = tpl.base_cost + 0.04 * (max(frag_sizes) / max(parent_heavy, 1))
                step = RouteStep(
                    product=canonical_smiles,
                    reactants=reactants,
                    reaction_smiles=f"{'.'.join(reactants)}>>{canonical_smiles}",
                    template_or_reference=f"{tpl.name}: {tpl.smarts}",
                    source="retro_star_smarts" if tpl.name not in RETRO_TEMPLATE_DELTAS else "template_fallback",
                )
                rank_key = (stock_misses, step_cost, max(frag_sizes))
                candidates.append((rank_key, step_cost, step))

        candidates.sort(key=lambda item: item[0])
        return [(c, st) for _, c, st in candidates[: self.max_branch_per_node]]

    def _compute_v0_heuristic(
        self,
        nm: NormalizedMol,
        v0_cache: dict[str, float],
    ) -> float:
        """Calibrated non-negative synthesis value estimate `V_0(m)` from `SCScore` & `SAscore`."""
        if nm.inchikey in v0_cache:
            return v0_cache[nm.inchikey]
        sa = compute_sascore(nm.canonical_smiles, mol=nm.mol)
        sc = compute_scscore(nm.canonical_smiles, weight_path=self.scscore_weights, mol=nm.mol)
        sa_term = max(0.0, ((sa or 2.5) - 1.8) * 0.22)
        if sc is not None:
            sc_term = max(0.0, (sc - 1.5) * 0.45)
            val = sc_term + sa_term
        else:
            val = max(0.0, ((sa or 2.5) - 1.7) * 0.35)
        v0_cache[nm.inchikey] = round(val, 6)
        return v0_cache[nm.inchikey]

    def search(
        self,
        target_smiles: str,
        target_inchikey: str,
        *,
        max_steps: int = 5,
        route_backend: str = "auto",
        offer_lookup: OfferLookupFn,
        lookup_backend: Optional[ReactionLookupBackend] = None,
        template_backend: Optional[TemplateFallbackBackend] = None,
        offer_cache: Optional[dict[str, list[Offer]]] = None,
        target_mol: Optional[Chem.Mol] = None,
        currency: str = "USD",
        desalt_fallback: bool = True,
        default_step_yield: float = 1.0,
        stoich_excess: float = 1.0,
        target_scale_g: float = 1.0,
        max_pack_multiple: Optional[float] = None,
        offer_mol_weights: Optional[Mapping[str, float]] = None,
    ) -> list[Route]:
        """Run `Retro*` best-first AND-OR search up to `max_steps` and optimize discrete MOQ outlay."""
        if max_steps <= 0:
            return []
        eff_max_steps = min(int(max_steps), 6)
        cache: dict[str, list[Offer]] = offer_cache if offer_cache is not None else {}
        nm_cache: dict[str, Optional[NormalizedMol]] = {}
        v0_cache: dict[str, float] = {}
        expansion_memo: dict[str, list[tuple[float, RouteStep]]] = {}

        def _get_nm(smi: str) -> Optional[NormalizedMol]:
            if smi not in nm_cache:
                nm_cache[smi] = normalize_smiles(smi, desalt=False)
            return nm_cache[smi]

        def _get_offers(nm: NormalizedMol) -> list[Offer]:
            return _resolve_leaf_offers(
                nm, cache, offer_lookup, desalt_fallback=desalt_fallback
            )

        def _is_in_stock(smi: str) -> bool:
            nm = _get_nm(smi)
            if nm is None or nm.inchikey == target_inchikey:
                return False
            return len(_get_offers(nm)) > 0

        def _leaf_price_tiebreak(nm: NormalizedMol) -> float:
            offs = _get_offers(nm)
            prices = [
                off.price / g
                for off in offs
                if off.price is not None
                and off.currency == currency
                and (g := offer_pack_grams(off)) is not None
                and g > 0
            ]
            if not prices:
                return 0.02
            return min(0.25, 0.01 * math.log1p(min(prices)))

        root_nm = _get_nm(target_smiles)
        if root_nm is None:
            return []
        root_mol = target_mol if target_mol is not None else root_nm.mol
        root_v0 = self._compute_v0_heuristic(root_nm, v0_cache)

        next_or_id = 0
        next_and_id = 0

        def _make_or_node(
            nm: NormalizedMol,
            depth: int,
            ancestors: frozenset[str],
            parent: Optional[_AndNode],
        ) -> _OrNode:
            nonlocal next_or_id
            nid = next_or_id
            next_or_id += 1
            stocked = depth > 0 and len(_get_offers(nm)) > 0
            if stocked:
                bonus = _leaf_price_tiebreak(nm)
                return _OrNode(
                    node_id=nid,
                    smiles=nm.canonical_smiles,
                    inchikey=nm.inchikey,
                    mol=nm.mol,
                    depth=depth,
                    ancestors=ancestors,
                    parent=parent,
                    in_stock=True,
                    v0=0.0,
                    leaf_cost_bonus=bonus,
                    expanded=True,
                    succ=True,
                    v_self=bonus,
                    min_steps=0,
                )
            v0 = root_v0 if depth == 0 else self._compute_v0_heuristic(nm, v0_cache)
            return _OrNode(
                node_id=nid,
                smiles=nm.canonical_smiles,
                inchikey=nm.inchikey,
                mol=root_mol if depth == 0 else nm.mol,
                depth=depth,
                ancestors=ancestors,
                parent=parent,
                in_stock=False,
                v0=v0,
                expanded=False,
                succ=False,
                v_self=v0,
                min_steps=1,
            )

        root = _make_or_node(
            root_nm,
            depth=0,
            ancestors=frozenset({target_inchikey}),
            parent=None,
        )

        use_lookup = (
            route_backend in {"auto", "lookup", "retro_star"}
            and lookup_backend is not None
            and lookup_backend.available
        )
        use_templates = route_backend in {"auto", "template", "retro_star"}
        tb = template_backend or TemplateFallbackBackend(max_branch_per_node=self.max_branch_per_node)

        def _candidate_disconnections(node: _OrNode) -> list[tuple[float, RouteStep]]:
            memo_key = f"{node.inchikey}:{'r' if node.depth == 0 else 'i'}"
            if memo_key in expansion_memo:
                return expansion_memo[memo_key]

            proposals: list[tuple[float, RouteStep]] = []
            seen_keys: set[tuple[str, ...]] = set()

            if use_lookup and lookup_backend is not None:
                l_steps = lookup_backend.one_step(
                    node.inchikey, node.smiles, limit=self.max_branch_per_node
                )
                for rank_idx, st in enumerate(l_steps):
                    key = tuple(sorted(st.reactants))
                    if not key or key in seen_keys or node.smiles in st.reactants:
                        continue
                    seen_keys.add(key)
                    # Calibrated negative log-prior cost favoring empirical patent/ORD reactions
                    rxn_cost = 0.72 + 0.06 * min(rank_idx, 8)
                    proposals.append((rxn_cost, st))

            if use_templates:
                if self.use_extended_templates:
                    t_steps = self.expand_extended_smarts(
                        node.smiles,
                        mol=node.mol,
                        is_stock_fn=_is_in_stock,
                        allow_unstocked_fgi=(node.depth == 0 and eff_max_steps == 1),
                    )
                    for step_cost, st in t_steps:
                        key = tuple(sorted(st.reactants))
                        if key in seen_keys:
                            continue
                        seen_keys.add(key)
                        proposals.append((step_cost, st))
                else:
                    for rank_idx, st in enumerate(tb.one_step(node.smiles, mol=node.mol)):
                        key = tuple(sorted(st.reactants))
                        if key in seen_keys:
                            continue
                        seen_keys.add(key)
                        proposals.append((1.15 + 0.05 * rank_idx, st))

            proposals = proposals[: self.max_branch_per_node]
            expansion_memo[memo_key] = proposals
            return proposals

        def _update_upward(start_or: _OrNode) -> None:
            curr: Optional[_OrNode] = start_or
            while curr is not None:
                if curr.in_stock:
                    curr.succ = True
                    curr.v_self = curr.leaf_cost_bonus
                    curr.min_steps = 0
                elif not curr.expanded:
                    curr.succ = False
                    curr.v_self = math.inf if curr.is_dead else curr.v0
                    curr.min_steps = 999 if curr.is_dead else 1
                else:
                    if not curr.children:
                        curr.is_dead = True
                        curr.succ = False
                        curr.v_self = math.inf
                        curr.min_steps = 999
                    else:
                        curr.succ = any(r.succ for r in curr.children)
                        curr.v_self = min(r.v_self for r in curr.children)
                        curr.min_steps = min(r.min_steps for r in curr.children)
                        if math.isinf(curr.v_self):
                            curr.is_dead = True

                and_parent = curr.parent
                if and_parent is None:
                    break
                and_parent.succ = all(c.succ for c in and_parent.children)
                if any(math.isinf(c.v_self) for c in and_parent.children):
                    and_parent.v_self = math.inf
                    and_parent.min_steps = 999
                else:
                    and_parent.v_self = and_parent.cost + sum(c.v_self for c in and_parent.children)
                    and_parent.min_steps = 1 + sum(c.min_steps for c in and_parent.children)
                curr = and_parent.parent

        def _select_best_open_node() -> Optional[_OrNode]:
            """Select open OR-node minimizing target-rooted `V_t(m | T)` subject to `<= eff_max_steps`."""
            best_cand: Optional[tuple[tuple[float, int, int, int], _OrNode]] = None

            def _walk(or_node: _OrNode, anc_v: float, anc_steps: int) -> None:
                nonlocal best_cand
                if or_node.in_stock or or_node.is_dead:
                    return
                total_min_steps = anc_steps + or_node.min_steps
                if total_min_steps > eff_max_steps:
                    return
                v_target = anc_v + or_node.v_self
                if not or_node.expanded:
                    if or_node.depth < eff_max_steps:
                        key = (v_target, total_min_steps, or_node.depth, or_node.node_id)
                        if best_cand is None or key < best_cand[0]:
                            best_cand = (key, or_node)
                    return

                for rxn_node in or_node.children:
                    if math.isinf(rxn_node.v_self):
                        continue
                    rxn_total_steps = anc_steps + rxn_node.min_steps
                    if rxn_total_steps > eff_max_steps:
                        continue
                    for child in rxn_node.children:
                        if child.in_stock or child.is_dead:
                            continue
                        sibling_v = rxn_node.v_self - child.v_self
                        sibling_steps = rxn_node.min_steps - child.min_steps
                        _walk(child, anc_v + sibling_v, anc_steps + sibling_steps)

            _walk(root, 0.0, 0)
            return best_cand[1] if best_cand is not None else None

        expansions = 0
        post_solve_used = 0
        while expansions < self.max_expansions:
            if root.succ:
                if post_solve_used >= self.post_solve_expansions:
                    break
                post_solve_used += 1

            chosen = _select_best_open_node()
            if chosen is None:
                break

            expansions += 1
            chosen.expanded = True
            cand_rxns = _candidate_disconnections(chosen)
            created_children: list[_AndNode] = []

            for rxn_cost, st in cand_rxns:
                child_nms: list[NormalizedMol] = []
                cycle_or_invalid = False
                for r_smi in st.reactants:
                    r_nm = _get_nm(r_smi)
                    if r_nm is None or r_nm.inchikey in chosen.ancestors:
                        cycle_or_invalid = True
                        break
                    child_nms.append(r_nm)
                if cycle_or_invalid or not child_nms:
                    continue

                and_node = _AndNode(
                    rxn_id=next_and_id,
                    step=st,
                    cost=rxn_cost,
                    parent=chosen,
                )
                next_and_id += 1

                next_anc = chosen.ancestors | {r_nm.inchikey for r_nm in child_nms}
                and_node.children = [
                    _make_or_node(
                        r_nm,
                        depth=chosen.depth + 1,
                        ancestors=chosen.ancestors | {r_nm.inchikey},
                        parent=and_node,
                    )
                    for r_nm in child_nms
                ]
                _ = next_anc
                and_node.succ = all(c.succ for c in and_node.children)
                and_node.v_self = and_node.cost + sum(c.v_self for c in and_node.children)
                and_node.min_steps = 1 + sum(c.min_steps for c in and_node.children)
                created_children.append(and_node)

            chosen.children = created_children
            _update_upward(chosen)

        return self._extract_routes(
            root,
            target_smiles=root_nm.canonical_smiles,
            target_inchikey=target_inchikey,
            max_steps=eff_max_steps,
            offer_cache=cache,
            offer_lookup=offer_lookup,
            currency=currency,
            desalt_fallback=desalt_fallback,
            default_step_yield=default_step_yield,
            stoich_excess=stoich_excess,
            target_scale_g=target_scale_g,
            max_pack_multiple=max_pack_multiple,
            offer_mol_weights=offer_mol_weights,
        )

    def _extract_routes(
        self,
        root: _OrNode,
        *,
        target_smiles: str,
        target_inchikey: str,
        max_steps: int,
        offer_cache: dict[str, list[Offer]],
        offer_lookup: OfferLookupFn,
        currency: str,
        desalt_fallback: bool,
        default_step_yield: float,
        stoich_excess: float,
        target_scale_g: float,
        max_pack_multiple: Optional[float],
        offer_mol_weights: Optional[Mapping[str, float]],
    ) -> list[Route]:
        """Extract top-K solved routes (or best partial disconnections if unsolved) and optimize MOQ."""

        def _enumerate_solved_subtrees(
            or_node: _OrNode, k_limit: int
        ) -> list[tuple[float, int, list[RouteStep], list[str]]]:
            if or_node.in_stock and or_node.depth > 0:
                return [(or_node.leaf_cost_bonus, 0, [], [or_node.smiles])]
            results: list[tuple[float, int, list[RouteStep], list[str]]] = []
            solved_and = sorted(
                (r for r in or_node.children if r.succ),
                key=lambda r: (r.min_steps, r.v_self, r.rxn_id),
            )
            for rxn in solved_and[:k_limit]:
                # Cartesian combine top branches per reactant child
                combos: list[tuple[float, int, list[RouteStep], list[str]]] = [
                    (rxn.cost, 0, [rxn.step], [])
                ]
                possible = True
                for child in rxn.children:
                    c_subs = _enumerate_solved_subtrees(child, max(2, k_limit))
                    if not c_subs:
                        possible = False
                        break
                    next_combos: list[tuple[float, int, list[RouteStep], list[str]]] = []
                    for acc_v, acc_d, acc_steps, acc_leaves in combos:
                        for cv, cd, csteps, cleaves in c_subs:
                            merged_steps = acc_steps + csteps
                            if len(merged_steps) <= max_steps:
                                next_combos.append(
                                    (
                                        acc_v + cv,
                                        max(acc_d, cd + 1),
                                        merged_steps,
                                        acc_leaves + cleaves,
                                    )
                                )
                    next_combos.sort(key=lambda t: (len(t[2]), t[0]))
                    combos = next_combos[:k_limit]
                    if not combos:
                        possible = False
                        break
                if possible:
                    results.extend(combos)
            results.sort(key=lambda t: (len(t[2]), t[0]))
            return results[:k_limit]

        candidates_raw: list[tuple[float, int, list[RouteStep], list[str]]] = []
        if root.succ:
            candidates_raw = _enumerate_solved_subtrees(root, self.top_k_routes * 2)

        # If unsolved within max_steps, emit best 1-step candidate disconnections so callers inspect partial leaves
        if not candidates_raw and root.children:
            partial_rxns = sorted(root.children, key=lambda r: (r.cost, len(r.step.reactants)))
            for rxn in partial_rxns[:3]:
                candidates_raw.append((rxn.cost, 1, [rxn.step], list(rxn.step.reactants)))

        built_routes: list[tuple[tuple, Route]] = []
        seen_route_sigs: set[tuple[tuple[str, tuple[str, ...]], ...]] = set()

        for v_val, tree_depth, steps, leaves in candidates_raw:
            if not steps or len(steps) > max_steps:
                continue
            sig = tuple((s.product, tuple(sorted(s.reactants))) for s in steps)
            if sig in seen_route_sigs:
                continue
            seen_route_sigs.add(sig)

            leaves_purchasable: list[bool] = []
            leaf_offers_by_ik: dict[str, list[Offer]] = {}
            for leaf_smi in leaves:
                lnm = normalize_smiles(leaf_smi)
                if lnm is None:
                    leaves_purchasable.append(False)
                    continue
                offs = _resolve_leaf_offers(
                    lnm, offer_cache, offer_lookup, desalt_fallback=desalt_fallback
                )
                leaf_offers_by_ik[lnm.inchikey] = offs
                leaves_purchasable.append(len(offs) > 0)

            breakdown = estimate_route_cost_bounds(
                target_smiles,
                leaves,
                leaf_offers_by_ik,
                currency=currency,
                steps=steps,
                default_step_yield=default_step_yield,
                stoich_excess=stoich_excess,
                target_scale_g=target_scale_g,
                max_pack_multiple=max_pack_multiple,
                offer_mol_weights=offer_mol_weights,
            )
            moq_opt = optimize_route_moq_outlay(
                target_smiles,
                leaves,
                leaf_offers_by_ik,
                steps=steps,
                currency=currency,
                default_step_yield=default_step_yield,
                stoich_excess=stoich_excess,
                target_scale_g=target_scale_g,
                offer_mol_weights=offer_mol_weights,
            )
            min_outlay = (
                moq_opt.pooled_outlay
                if moq_opt.pooled_outlay is not None
                else breakdown.min_outlay_cost
            )

            sources_used = ", ".join(dict.fromkeys(s.source for s in steps))
            score = round(1.0 / (1.0 + max(0.0, float(v_val))), 4)
            route = Route(
                target_inchikey=target_inchikey,
                n_steps=len(steps),
                tree_depth=tree_depth,
                steps=steps,
                leaves=leaves,
                leaves_purchasable=leaves_purchasable,
                method=f"retro_star (best-first AND-OR search: {sources_used})",
                score=score,
                est_cost=breakdown.est_cost,
                est_cost_method=route_cost_method(
                    breakdown.skipped,
                    default_step_yield=default_step_yield,
                    stoich_excess=stoich_excess,
                ),
                est_cost_complete=breakdown.est_cost_complete,
                est_cost_yield_adjusted=breakdown.yield_adjusted_cost,
                min_outlay_cost=min_outlay,
                step_yield_assumed=default_step_yield,
            )

            all_ok = 0 if (leaves and all(leaves_purchasable)) else 1
            cost_tier = (
                0
                if (route.est_cost_complete and route.est_cost is not None)
                else (1 if route.est_cost is not None else 2)
            )
            outlay_val = (
                float(min_outlay.value)
                if (min_outlay is not None and min_outlay.value is not None and moq_opt.complete)
                else math.inf
            )
            est_val = (
                float(route.est_cost.value)
                if (route.est_cost is not None and route.est_cost.value is not None)
                else math.inf
            )
            sort_key = (all_ok, route.n_steps, tree_depth, cost_tier, outlay_val, est_val, -score)
            built_routes.append((sort_key, route))

        built_routes.sort(key=lambda item: item[0])
        return [r for _, r in built_routes[: self.top_k_routes]]


__all__ = [
    "EXTENDED_RETRO_STAR_TEMPLATES",
    "MOQOptimizationResult",
    "PackAllocation",
    "RETRO_STAR_TEMPLATE_DELTAS",
    "RetroStarPlanner",
    "UNIMOLECULAR_RETRO_STAR_TEMPLATES",
    "optimize_discrete_pack_order",
    "optimize_route_moq_outlay",
]
