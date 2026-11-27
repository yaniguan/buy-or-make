"""Retrosynthetic route backends: reaction lookup (USPTO/ORD), template fallback, optional AiZynthFinder."""

from __future__ import annotations

import json
import os
import sqlite3
import threading
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Optional

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .cache import open_readonly_sqlite
from .canonical import NormalizedMol, normalize_smiles
from .cost import EST_COST_METHOD, estimate_route_cost_bounds, route_cost_method
from .schema import Offer, Route, RouteStep

DEFAULT_REACTION_DB = Path("data/reactions.sqlite")

# Curated open SMARTS retrosynthesis disconnections (product >> reactants).
# Every bimolecular template enforces `!@` on the cleaved bond so endocyclic bonds are never
# split into atom-duplicating fragments; monocyclic lactones/lactams use dedicated single-component
# ring-opening SMARTS (`@` with component grouping). Williamson SN2 templates restrict electrophile
# carbons to methyl, primary, or unhindered secondary sp3 centers (`H3,H2,H1`, non-tertiary).
RETRO_TEMPLATES: list[tuple[str, str]] = [
    (
        "amide_coupling_secondary",
        "[#6:1][C:2](=[O:3])!@[NH1:4][#6:5]>>[#6:1][C:2](=[O:3])O.[NH2:4][#6:5]",
    ),
    (
        "amide_coupling_tertiary",
        "[#6:1][C:2](=[O:3])!@[N:4]([#6:5])[#6:6]>>[#6:1][C:2](=[O:3])O.[NH1:4]([#6:5])[#6:6]",
    ),
    (
        "esterification",
        "[#6:1][C:2](=[O:3])!@[O:4][#6:5]>>[#6:1][C:2](=[O:3])O.[OH1:4][#6:5]",
    ),
    (
        "sulfonamide_secondary",
        "[#6:1][S:2](=[O:3])(=[O:4])!@[NH1:5][#6:6]>>[#6:1][S:2](=[O:3])(=[O:4])Cl.[NH2:5][#6:6]",
    ),
    (
        "sulfonamide_tertiary",
        "[#6:1][S:2](=[O:3])(=[O:4])!@[N:5]([#6:6])[#6:7]>>[#6:1][S:2](=[O:3])(=[O:4])Cl.[NH1:5]([#6:6])[#6:7]",
    ),
    (
        "urea_formation",
        "[#6:1][NH1:2][C:3](=[O:4])!@[NH1:5][#6:6]>>[#6:1][N:2]=[C:3]=[O:4].[NH2:5][#6:6]",
    ),
    (
        "carbamate_formation",
        "[#6:1][O:2]!@[C:3](=[O:4])[NH1:5][#6:6]>>[#6:1][OH1:2].[#6:6][N:5]=[C:3]=[O:4]",
    ),
    (
        "suzuki_biaryl",
        "[c:1]-!@[c:2]>>[c:1]Br.[c:2]B(O)O",
    ),
    (
        "suzuki_aryl_alkenyl",
        "[c:1]-!@[CX3:2]=[CX3:3]>>[c:1]Br.[CX3:2](B(O)O)=[CX3:3]",
    ),
    (
        "williamson_aryl_ether",
        "[c:1][O:2]!@[CX4;H3,H2,H1;!$(C([#6])([#6])[#6]):3]>>[c:1][OH1:2].[CX4:3]Br",
    ),
    (
        "williamson_alkyl_ether",
        "[CX4:1][O:2]!@[CX4;H3,H2,H1;!$(C([#6])([#6])[#6]):3]>>[CX4:1][OH1:2].[CX4:3]Br",
    ),
    (
        "reductive_amination_secondary",
        "[#6:1][CH2:2]!@[NH1:3][#6:4]>>[#6:1][CH1:2]=O.[NH2:3][#6:4]",
    ),
    (
        "reductive_amination_tertiary",
        "[#6:1][CH2:2]!@[N:3]([#6:4])[#6:5]>>[#6:1][CH1:2]=O.[NH1:3]([#6:4])[#6:5]",
    ),
    (
        "buchwald_hartwig_cn",
        "[c:1]!@[NH1:2][c:3]>>[c:1]Br.[NH2:2][c:3]",
    ),
    (
        "buchwald_hartwig_cn_tertiary",
        "[c:1]!@[N:2]([#6:3])[#6:4]>>[c:1]Br.[NH1:2]([#6:3])[#6:4]",
    ),
    (
        "lactone_ring_opening",
        "([#6:1][C:2](=[O:3])@[O:4][#6:5])>>([#6:1][C:2](=[O:3])O.[OH1:4][#6:5])",
    ),
    (
        "lactam_ring_opening",
        "([#6:1][C:2](=[O:3])@[NH1:4][#6:5])>>([#6:1][C:2](=[O:3])O.[NH2:4][#6:5])",
    ),
]

UNIMOLECULAR_RING_OPENING_TEMPLATES: frozenset[str] = frozenset(
    {"lactone_ring_opening", "lactam_ring_opening"}
)

# Exact heavy-atom leaving-group delta (sum(reactant heavy atoms) - product heavy atoms) per template.
RETRO_TEMPLATE_DELTAS: dict[str, dict[str, int]] = {
    "amide_coupling_secondary": {"O": 1},
    "amide_coupling_tertiary": {"O": 1},
    "esterification": {"O": 1},
    "sulfonamide_secondary": {"Cl": 1},
    "sulfonamide_tertiary": {"Cl": 1},
    "urea_formation": {},
    "carbamate_formation": {},
    "suzuki_biaryl": {"Br": 1, "B": 1, "O": 2},
    "suzuki_aryl_alkenyl": {"Br": 1, "B": 1, "O": 2},
    "williamson_aryl_ether": {"Br": 1},
    "williamson_alkyl_ether": {"Br": 1},
    "reductive_amination_secondary": {"O": 1},
    "reductive_amination_tertiary": {"O": 1},
    "buchwald_hartwig_cn": {"Br": 1},
    "buchwald_hartwig_cn_tertiary": {"Br": 1},
    "lactone_ring_opening": {"O": 1},
    "lactam_ring_opening": {"O": 1},
}

_CATALYST_METALS = {"Pd", "Pt", "Rh", "Ir", "Ru", "Ni", "Au"}
_INORGANIC_COUNTER_METALS = {"Li", "Na", "K", "Rb", "Cs", "Ag"}

# Pre-compiled chemoselectivity / homopolymerization SMARTS guards
_PAT_UNPROTECTED_AMINE = Chem.MolFromSmarts(
    "[N;!a;H2,H1;!$(NC(=[O,S]));!$(NS(=O)=O);!$(N=[C,N,O]);!$(N[N,O])]"
)
_PAT_UNPROTECTED_OH = Chem.MolFromSmarts(
    "[OH1;!$(OC(=[O,S]));!$(OS(=O));!$(OB);!$(OP);!$(ON)]"
)
_PAT_CARBOXYLIC_ACID = Chem.MolFromSmarts("[CX3](=O)[OH1]")
_PAT_REACTIVE_ARYL_VINYL_HAL = Chem.MolFromSmarts("[c,CX3][Br,I]")
_PAT_BORONIC_ACID = Chem.MolFromSmarts("[#6]B([OH])[OH]")
_PAT_ALDEHYDE = Chem.MolFromSmarts("[CX3H1](=O)[#6,#1]")


def _has_chemoselectivity_conflict(tpl_name: str, frag_mols: list[Chem.Mol]) -> bool:
    """Return True if disconnected reactants carry unprotected self-reacting or bis-coupling groups."""
    if len(frag_mols) < 2:
        return False
    f0, f1 = frag_mols[0], frag_mols[1]

    if tpl_name in {"amide_coupling_secondary", "amide_coupling_tertiary"}:
        # Acyl donor (f0) must not carry an unprotected primary/secondary amine (AB self-condensation)
        # or multiple carboxylic acids; amine nucleophile (f1) must not have >=2 competing free amines.
        if f0.HasSubstructMatch(_PAT_UNPROTECTED_AMINE):
            return True
        if len(f0.GetSubstructMatches(_PAT_CARBOXYLIC_ACID)) >= 2:
            return True
        if len(f1.GetSubstructMatches(_PAT_UNPROTECTED_AMINE)) >= 2:
            return True
        return False

    if tpl_name == "esterification":
        # Carboxylic acid donor (f0) must not carry unprotected amine or alcohol/phenol (lactonization/polyester)
        if f0.HasSubstructMatch(_PAT_UNPROTECTED_AMINE) or f0.HasSubstructMatch(_PAT_UNPROTECTED_OH):
            return True
        return False

    if tpl_name in {"sulfonamide_secondary", "sulfonamide_tertiary"}:
        if f0.HasSubstructMatch(_PAT_UNPROTECTED_AMINE) or f0.HasSubstructMatch(_PAT_UNPROTECTED_OH):
            return True
        if len(f1.GetSubstructMatches(_PAT_UNPROTECTED_AMINE)) >= 2:
            return True
        return False

    if tpl_name == "urea_formation":
        # Isocyanate (f0) self-reacts with unprotected amines or alcohols
        if f0.HasSubstructMatch(_PAT_UNPROTECTED_AMINE) or f0.HasSubstructMatch(_PAT_UNPROTECTED_OH):
            return True
        return False

    if tpl_name == "carbamate_formation":
        # Isocyanate (f1) self-reacts with unprotected amines or alcohols
        if f1.HasSubstructMatch(_PAT_UNPROTECTED_AMINE) or f1.HasSubstructMatch(_PAT_UNPROTECTED_OH):
            return True
        return False

    if tpl_name in {"suzuki_biaryl", "suzuki_aryl_alkenyl"}:
        # Electrophile (f0 = Ar-Br) must not carry >=2 reactive [c,CX3][Br,I] sites (bis-coupling) or -B(O)O;
        # Boronic acid (f1) must not carry reactive [c,CX3][Br,I] (AB Suzuki homopolymerization) or >=2 -B(O)O.
        if len(f0.GetSubstructMatches(_PAT_REACTIVE_ARYL_VINYL_HAL)) >= 2:
            return True
        if f0.HasSubstructMatch(_PAT_BORONIC_ACID):
            return True
        if f1.HasSubstructMatch(_PAT_REACTIVE_ARYL_VINYL_HAL):
            return True
        if len(f1.GetSubstructMatches(_PAT_BORONIC_ACID)) >= 2:
            return True
        return False

    if tpl_name in {"buchwald_hartwig_cn", "buchwald_hartwig_cn_tertiary"}:
        # Aryl bromide (f0) must not carry unprotected amine (AB self-amination) or >=2 reactive Ar-Br/I sites
        if f0.HasSubstructMatch(_PAT_UNPROTECTED_AMINE):
            return True
        if len(f0.GetSubstructMatches(_PAT_REACTIVE_ARYL_VINYL_HAL)) >= 2:
            return True
        return False

    if tpl_name in {"reductive_amination_secondary", "reductive_amination_tertiary"}:
        # Aldehyde (f0) must not carry an unprotected amine (self-imine polymerization) or >=2 aldehydes
        if f0.HasSubstructMatch(_PAT_UNPROTECTED_AMINE):
            return True
        if len(f0.GetSubstructMatches(_PAT_ALDEHYDE)) >= 2:
            return True
        return False

    if tpl_name in {"williamson_aryl_ether", "williamson_alkyl_ether"}:
        # Alkyl bromide (f1) must not be tertiary/bridgehead (E2 elimination) or carry self-alkylating -NH/-OH
        if f1.HasSubstructMatch(_PAT_UNPROTECTED_AMINE) or f1.HasSubstructMatch(_PAT_UNPROTECTED_OH):
            return True
        ring_info = f1.GetRingInfo()
        for atom in f1.GetAtoms():
            if atom.GetAtomicNum() == 35:  # Br
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6:
                        if nbr.GetTotalNumHs() == 0:
                            return True
                        if ring_info.NumAtomRings(nbr.GetIdx()) >= 2:
                            return True
        return False

    return False


def _heavy_element_counts(mol: Chem.Mol) -> Counter[str]:
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def parse_reaction_smiles_for_indexing(
    rxn_smi: str,
) -> list[tuple[NormalizedMol, list[str], list[str]]]:
    """Extract `(product_mol, stoichiometric_reactants, spectator_agents)` per product in `rxn_smi`."""
    if not rxn_smi or ">>" not in rxn_smi:
        return []
    lhs, rhs = rxn_smi.split(">>", 1)
    lhs_parts = lhs.split(">")
    raw_reactants = [r.strip() for r in lhs_parts[0].split(".") if r.strip()]
    raw_explicit_agents = [
        a.strip() for part in lhs_parts[1:] for a in part.split(".") if a.strip()
    ]
    raw_products = [p.strip().split()[0] for p in rhs.split(".") if p.strip()]
    if not raw_reactants or not raw_products:
        return []

    parsed_reactants: list[tuple[Chem.Mol, set[int], NormalizedMol]] = []
    for r_raw in raw_reactants:
        r_mol = Chem.MolFromSmiles(r_raw)
        if r_mol is None:
            continue
        r_maps = {
            a.GetAtomMapNum()
            for a in r_mol.GetAtoms()
            if a.GetAtomicNum() > 1 and a.GetAtomMapNum() > 0
        }
        r_nm = normalize_smiles(r_raw)
        if r_nm is None:
            continue
        parsed_reactants.append((r_mol, r_maps, r_nm))

    if not parsed_reactants:
        return []

    explicit_agent_smis = [
        nm.canonical_smiles for a in raw_explicit_agents if (nm := normalize_smiles(a)) is not None
    ]

    extracted: list[tuple[NormalizedMol, list[str], list[str]]] = []
    for p_raw in raw_products:
        p_mol = Chem.MolFromSmiles(p_raw)
        if p_mol is None:
            continue
        prod_maps = {
            a.GetAtomMapNum()
            for a in p_mol.GetAtoms()
            if a.GetAtomicNum() > 1 and a.GetAtomMapNum() > 0
        }
        p_nm = normalize_smiles(p_raw)
        if p_nm is None:
            continue
        prod_elements = {a.GetSymbol() for a in p_mol.GetAtoms() if a.GetAtomicNum() > 1}

        stoich: list[str] = []
        agents: list[str] = list(explicit_agent_smis)

        if prod_maps:
            for _, r_maps, r_nm in parsed_reactants:
                if r_maps & prod_maps:
                    stoich.append(r_nm.canonical_smiles)
                else:
                    agents.append(r_nm.canonical_smiles)
        else:
            for r_mol, _, r_nm in parsed_reactants:
                r_elems = {a.GetSymbol() for a in r_mol.GetAtoms() if a.GetAtomicNum() > 1}
                has_foreign_catalyst_metal = bool((r_elems & _CATALYST_METALS) - prod_elements)
                n_carbons = sum(1 for a in r_mol.GetAtoms() if a.GetAtomicNum() == 6)
                is_inorganic_spectator = (
                    n_carbons == 0 and bool((r_elems & _INORGANIC_COUNTER_METALS) - prod_elements)
                )
                if has_foreign_catalyst_metal or is_inorganic_spectator:
                    agents.append(r_nm.canonical_smiles)
                else:
                    stoich.append(r_nm.canonical_smiles)

        if stoich:
            extracted.append((p_nm, stoich, agents))
    return extracted


@dataclass(frozen=True)
class _CompiledTemplate:
    name: str
    smarts: str
    rxn: rdChemReactions.ChemicalReaction
    heavy_delta: dict[str, int]


_COMPILED_TEMPLATES: Optional[list[_CompiledTemplate]] = None


def _get_compiled_templates() -> list[_CompiledTemplate]:
    global _COMPILED_TEMPLATES
    if _COMPILED_TEMPLATES is None:
        _COMPILED_TEMPLATES = [
            _CompiledTemplate(
                name=n,
                smarts=s,
                rxn=rdChemReactions.ReactionFromSmarts(s),
                heavy_delta=dict(RETRO_TEMPLATE_DELTAS.get(n, {})),
            )
            for n, s in RETRO_TEMPLATES
        ]
    return _COMPILED_TEMPLATES


class ReactionLookupBackend:
    """Exact product-InChIKey lookup in locally indexed USPTO (CC0) / ORD (CC-BY-SA-4.0) reactions."""

    def __init__(self, db_path: Path | str = DEFAULT_REACTION_DB) -> None:
        self.db_path = Path(db_path)
        self._lock = threading.RLock()
        self._conn: Optional[sqlite3.Connection] = None
        if self.db_path.exists():
            self._conn = open_readonly_sqlite(self.db_path)

    @property
    def available(self) -> bool:
        with self._lock:
            return self._conn is not None

    def one_step(self, inchikey: str, canonical_smiles: str, limit: int = 10) -> list[RouteStep]:
        with self._lock:
            if self._conn is None:
                return []
            rows = self._conn.execute(
                "SELECT reactants_json, reaction_smiles, reference, source "
                "FROM reactions WHERE product_inchikey = ? LIMIT ?",
                (inchikey, limit),
            ).fetchall()
        steps: list[RouteStep] = []
        seen: set[tuple[str, ...]] = set()
        for r_json, rxn_smi, ref, src in rows:
            reactants = json.loads(r_json)
            norm_reactants: list[str] = []
            valid = True
            for r in reactants:
                nm = normalize_smiles(r)
                if nm is None:
                    valid = False
                    break
                norm_reactants.append(nm.canonical_smiles)
            if not valid or not norm_reactants:
                continue
            key = tuple(sorted(norm_reactants))
            if key in seen:
                continue
            seen.add(key)
            steps.append(
                RouteStep(
                    product=canonical_smiles,
                    reactants=norm_reactants,
                    reaction_smiles=rxn_smi,
                    template_or_reference=ref,
                    source=src,
                )
            )
        return steps

    def close(self) -> None:
        with self._lock:
            if self._conn is not None:
                self._conn.close()
                self._conn = None


class TemplateFallbackBackend:
    """Built-in 1- and 2-step heuristic SMARTS retrosynthetic disconnection engine."""

    def __init__(self, max_branch_per_node: int = 12) -> None:
        self.max_branch_per_node = max_branch_per_node
        self.templates = _get_compiled_templates()

    def one_step(
        self,
        canonical_smiles: str,
        *,
        mol: Optional[Chem.Mol] = None,
    ) -> list[RouteStep]:
        parsed = mol if mol is not None else Chem.MolFromSmiles(canonical_smiles)
        if parsed is None:
            return []
        parent_heavy = parsed.GetNumHeavyAtoms()
        parent_counts = _heavy_element_counts(parsed)
        candidates: list[tuple[int, RouteStep]] = []
        seen_sets: set[tuple[str, ...]] = set()

        for tpl in self.templates:
            try:
                outcomes = tpl.rxn.RunReactants((parsed,))
            except Exception:
                continue
            expected_counts = Counter(parent_counts)
            for sym, delta in tpl.heavy_delta.items():
                expected_counts[sym] += delta
                if expected_counts[sym] <= 0:
                    del expected_counts[sym]

            for prod_tuple in outcomes:
                reactants: list[str] = []
                frag_mols: list[Chem.Mol] = []
                frag_heavy_sizes: list[int] = []
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
                    frag_heavy_sizes.append(nm.num_heavy_atoms)
                    frag_counts.update(_heavy_element_counts(nm.mol))
                    reactants.append(nm.canonical_smiles)

                if not ok or not reactants:
                    continue
                if frag_counts != expected_counts:
                    continue
                if any(sz >= parent_heavy + 3 for sz in frag_heavy_sizes):
                    continue
                if _has_chemoselectivity_conflict(tpl.name, frag_mols):
                    continue
                key = tuple(sorted(reactants))
                if key in seen_sets:
                    continue
                seen_sets.add(key)
                rxn_smi = f"{'.'.join(reactants)}>>{canonical_smiles}"
                step = RouteStep(
                    product=canonical_smiles,
                    reactants=reactants,
                    reaction_smiles=rxn_smi,
                    template_or_reference=f"{tpl.name}: {tpl.smarts}",
                    source="template_fallback",
                )
                max_frag = max(frag_heavy_sizes)
                candidates.append((max_frag, step))

        candidates.sort(key=lambda item: (item[0], len(item[1].reactants)))
        return [step for _, step in candidates[: self.max_branch_per_node]]


def _select_if_configured(policy_or_stock: object) -> None:
    """Call `select_all()` only when the policy/stock container actually has configured items."""
    items = getattr(policy_or_stock, "items", None)
    if items is not None and len(items) == 0:
        return
    select_all = getattr(policy_or_stock, "select_all", None)
    if callable(select_all):
        select_all()


class AiZynthFinderBackend:
    """Optional AiZynthFinder wrapper (MIT code, CC-BY-4.0 policy models)."""

    def __init__(self, config_path: Optional[Path | str] = None) -> None:
        cfg = config_path or os.environ.get("AIZYNTH_CONFIG")
        self.config_path = Path(cfg) if cfg else None
        self._finder = None

    @staticmethod
    def is_installed() -> bool:
        try:
            import aizynthfinder  # noqa: F401

            return True
        except ImportError:
            return False

    def is_ready(self) -> bool:
        return self.is_installed() and self.config_path is not None and self.config_path.exists()

    def search(self, target_smiles: str, target_inchikey: str, max_steps: int) -> list[Route]:
        if not self.is_ready():
            return []
        from aizynthfinder.aizynthfinder import AiZynthFinder  # type: ignore

        if self._finder is None:
            self._finder = AiZynthFinder(configfile=str(self.config_path))
            _select_if_configured(self._finder.expansion_policy)
            _select_if_configured(self._finder.filter_policy)
            _select_if_configured(self._finder.stock)

        self._finder.target_smiles = target_smiles
        self._finder.config.search.max_transforms = max_steps
        self._finder.tree_search()
        self._finder.build_routes()

        routes: list[Route] = []
        for route_dict in self._finder.routes[:5]:
            tree = route_dict.get("dict") or {}
            steps: list[RouteStep] = []
            leaves: list[str] = []
            leaves_ok: list[bool] = []
            _walk_aizynth_tree(tree, steps, leaves, leaves_ok)
            if not steps:
                continue
            raw_score = route_dict.get("score")
            if isinstance(raw_score, dict):
                raw_score = next(
                    (float(v) for v in raw_score.values() if isinstance(v, (int, float))),
                    None,
                )
            score_val = float(raw_score) if isinstance(raw_score, (int, float)) else None
            routes.append(
                Route(
                    target_inchikey=target_inchikey,
                    n_steps=len(steps),
                    steps=steps,
                    leaves=leaves,
                    leaves_purchasable=leaves_ok,
                    method="aizynthfinder (MCTS retrosynthesis)",
                    score=score_val,
                    est_cost=None,
                    est_cost_method=EST_COST_METHOD,
                    est_cost_complete=False,
                )
            )
        return routes


def _walk_aizynth_tree(
    node: dict,
    steps: list[RouteStep],
    leaves: list[str],
    leaves_ok: list[bool],
) -> None:
    """Walk an AiZynthFinder route tree, selecting the best reaction branch at multi-OR mol nodes."""
    node_type = node.get("type")
    if node_type == "mol":
        children = node.get("children") or []
        rxn_children = [c for c in children if isinstance(c, dict) and c.get("type") == "reaction"]
        if not rxn_children and children and isinstance(children[0], dict):
            rxn_children = [children[0]]
        if not rxn_children:
            leaves.append(node["smiles"])
            leaves_ok.append(bool(node.get("in_stock", False)))
            return

        best_branch: Optional[tuple[tuple[int, int, int], list[RouteStep], list[str], list[bool]]] = None
        for rxn_node in rxn_children:
            b_steps: list[RouteStep] = []
            b_leaves: list[str] = []
            b_ok: list[bool] = []
            reactant_nodes = [
                c for c in (rxn_node.get("children") or []) if isinstance(c, dict) and c.get("type") == "mol"
            ]
            if not reactant_nodes:
                continue
            reactant_smiles = [c["smiles"] for c in reactant_nodes]
            b_steps.append(
                RouteStep(
                    product=node["smiles"],
                    reactants=reactant_smiles,
                    reaction_smiles=f"{'.'.join(reactant_smiles)}>>{node['smiles']}",
                    template_or_reference=str(
                        (rxn_node.get("metadata") or {}).get("template_code", "aizynthfinder")
                    ),
                    source="aizynthfinder",
                )
            )
            for child in reactant_nodes:
                _walk_aizynth_tree(child, b_steps, b_leaves, b_ok)
            rank_key = (
                0 if (b_leaves and all(b_ok)) else 1,
                len(b_steps),
                len(b_leaves),
            )
            if best_branch is None or rank_key < best_branch[0]:
                best_branch = (rank_key, b_steps, b_leaves, b_ok)

        if best_branch is not None:
            steps.extend(best_branch[1])
            leaves.extend(best_branch[2])
            leaves_ok.extend(best_branch[3])
        else:
            leaves.append(node["smiles"])
            leaves_ok.append(bool(node.get("in_stock", False)))
    elif node_type == "reaction":
        reactant_nodes = [
            c for c in (node.get("children") or []) if isinstance(c, dict) and c.get("type") == "mol"
        ]
        for child in reactant_nodes:
            _walk_aizynth_tree(child, steps, leaves, leaves_ok)


OfferLookupFn = Callable[[str, str], list[Offer]]


def _resolve_leaf_offers(
    nm: NormalizedMol,
    offer_cache: dict[str, list[Offer]],
    offer_lookup: OfferLookupFn,
    *,
    desalt_fallback: bool = True,
) -> list[Offer]:
    if nm.inchikey not in offer_cache:
        offs = list(offer_lookup(nm.inchikey, nm.canonical_smiles))
        if not offs and desalt_fallback and (
            "." in nm.canonical_smiles
            or "+" in nm.canonical_smiles
            or "-" in nm.canonical_smiles
        ):
            p_nm = normalize_smiles(nm.canonical_smiles, desalt=True)
            if p_nm is not None and p_nm.inchikey != nm.inchikey:
                if p_nm.inchikey not in offer_cache:
                    offer_cache[p_nm.inchikey] = list(
                        offer_lookup(p_nm.inchikey, p_nm.canonical_smiles)
                    )
                offs = list(offer_cache[p_nm.inchikey])
        offer_cache[nm.inchikey] = offs
    return offer_cache[nm.inchikey]


def _score_route_for_selection(route: Route) -> tuple:
    all_ok = 0 if (route.leaves and all(route.leaves_purchasable)) else 1
    cost_tier = 2
    cost_val = 0.0
    yield_cost_val = 0.0
    if route.est_cost is not None and route.est_cost.value is not None:
        cost_tier = 0 if route.est_cost_complete else 1
        cost_val = route.est_cost.value
    if route.est_cost_yield_adjusted is not None and route.est_cost_yield_adjusted.value is not None:
        yield_cost_val = route.est_cost_yield_adjusted.value
    neg_score = -float(route.score) if route.score is not None else 0.0
    return (all_ok, route.n_steps, cost_tier, cost_val, yield_cost_val, neg_score, len(route.leaves))


def _build_route(
    target_smiles: str,
    target_inchikey: str,
    steps: list[RouteStep],
    leaves: list[str],
    method: str,
    offer_cache: dict[str, list[Offer]],
    offer_lookup: OfferLookupFn,
    *,
    score: Optional[float] = None,
    currency: str = "USD",
    desalt_fallback: bool = True,
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
    target_scale_g: float = 1.0,
    max_pack_multiple: Optional[float] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
) -> Route:
    leaves_purchasable: list[bool] = []
    leaf_offers_by_ik: dict[str, list[Offer]] = {}
    for leaf_smi in leaves:
        nm = normalize_smiles(leaf_smi)
        if nm is None:
            leaves_purchasable.append(False)
            continue
        offs = _resolve_leaf_offers(
            nm, offer_cache, offer_lookup, desalt_fallback=desalt_fallback
        )
        leaf_offers_by_ik[nm.inchikey] = offs
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
    return Route(
        target_inchikey=target_inchikey,
        n_steps=len(steps),
        steps=steps,
        leaves=leaves,
        leaves_purchasable=leaves_purchasable,
        method=method,
        score=score,
        est_cost=breakdown.est_cost,
        est_cost_method=route_cost_method(
            breakdown.skipped,
            default_step_yield=default_step_yield,
            stoich_excess=stoich_excess,
        ),
        est_cost_complete=breakdown.est_cost_complete,
        est_cost_yield_adjusted=breakdown.yield_adjusted_cost,
        min_outlay_cost=breakdown.min_outlay_cost,
        step_yield_assumed=default_step_yield,
    )


def _expand_lookup_routes(
    target_smiles: str,
    target_inchikey: str,
    max_steps: int,
    lookup_backend: ReactionLookupBackend,
    offer_cache: dict[str, list[Offer]],
    offer_lookup: OfferLookupFn,
    currency: str,
    *,
    desalt_fallback: bool = True,
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
    target_scale_g: float = 1.0,
    max_pack_multiple: Optional[float] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
) -> list[Route]:
    """Bounded multi-step reaction lookup expansion up to `min(max_steps, 3)` steps."""
    found: list[Route] = []
    step1_list = lookup_backend.one_step(target_inchikey, target_smiles)
    frontier: list[tuple[list[RouteStep], list[str], set[str]]] = []
    solved_at_current_depth = False

    for st1 in step1_list:
        sources_str = st1.source
        r1 = _build_route(
            target_smiles,
            target_inchikey,
            [st1],
            list(st1.reactants),
            f"reaction_lookup ({sources_str})",
            offer_cache,
            offer_lookup,
            currency=currency,
            desalt_fallback=desalt_fallback,
            default_step_yield=default_step_yield,
            stoich_excess=stoich_excess,
            target_scale_g=target_scale_g,
            max_pack_multiple=max_pack_multiple,
            offer_mol_weights=offer_mol_weights,
        )
        found.append(r1)
        if r1.leaves and all(r1.leaves_purchasable):
            solved_at_current_depth = True
        else:
            frontier.append(([st1], list(st1.reactants), {target_inchikey}))

    effective_depth = min(max_steps, 3)
    depth = 1
    while not solved_at_current_depth and depth < effective_depth and frontier:
        next_frontier: list[tuple[list[RouteStep], list[str], set[str]]] = []
        for steps_so_far, current_leaves, visited_iks in frontier[:6]:
            missing_idx = [
                idx
                for idx, leaf in enumerate(current_leaves)
                if not _is_purchasable(
                    leaf, offer_cache, offer_lookup, desalt_fallback=desalt_fallback
                )
            ]
            if len(missing_idx) != 1:
                continue
            m_idx = missing_idx[0]
            missing_leaf = current_leaves[m_idx]
            leaf_nm = normalize_smiles(missing_leaf)
            if leaf_nm is None or leaf_nm.inchikey in visited_iks:
                continue
            other_leaves = [r for i, r in enumerate(current_leaves) if i != m_idx]
            for st_next in lookup_backend.one_step(
                leaf_nm.inchikey, leaf_nm.canonical_smiles, limit=6
            ):
                new_steps = [*steps_so_far, st_next]
                new_leaves = [*other_leaves, *st_next.reactants]
                src_labels = ", ".join(dict.fromkeys(s.source for s in new_steps))
                r_multi = _build_route(
                    target_smiles,
                    target_inchikey,
                    new_steps,
                    new_leaves,
                    f"reaction_lookup ({src_labels})",
                    offer_cache,
                    offer_lookup,
                    currency=currency,
                    desalt_fallback=desalt_fallback,
                    default_step_yield=default_step_yield,
                    stoich_excess=stoich_excess,
                    target_scale_g=target_scale_g,
                    max_pack_multiple=max_pack_multiple,
                    offer_mol_weights=offer_mol_weights,
                )
                found.append(r_multi)
                if r_multi.leaves and all(r_multi.leaves_purchasable):
                    solved_at_current_depth = True
                    break
                next_frontier.append(
                    (new_steps, new_leaves, visited_iks | {leaf_nm.inchikey})
                )
        frontier = next_frontier
        depth += 1

    return found


def find_routes(
    target_smiles: str,
    target_inchikey: str,
    *,
    max_steps: int = 3,
    route_backend: str = "auto",
    planner: str = "lookup_smarts",
    offer_lookup: OfferLookupFn,
    lookup_backend: Optional[ReactionLookupBackend] = None,
    template_backend: Optional[TemplateFallbackBackend] = None,
    aizynth_backend: Optional[AiZynthFinderBackend] = None,
    offer_cache: Optional[dict[str, list[Offer]]] = None,
    target_mol: Optional[Chem.Mol] = None,
    currency: str = "USD",
    desalt_fallback: bool = True,
    default_step_yield: float = 1.0,
    stoich_excess: float = 1.0,
    target_scale_g: float = 1.0,
    max_pack_multiple: Optional[float] = None,
    offer_mol_weights: Optional[Mapping[str, float]] = None,
    scscore_weights: Optional[Path | str] = None,
    retro_star_max_expansions: int = 90,
) -> list[Route]:
    """Search for retrosynthetic routes up to max_steps using the selected backend and planner."""
    if max_steps <= 0:
        return []
    if planner not in {"lookup_smarts", "retro_star"}:
        raise ValueError(f"Unsupported planner {planner!r}; expected 'lookup_smarts' or 'retro_star'")

    cache = offer_cache if offer_cache is not None else {}

    if planner == "retro_star" or route_backend == "retro_star":
        from .retro_star import RetroStarPlanner
        from .scores import DEFAULT_SCSCORE_WEIGHTS

        rs_planner = RetroStarPlanner(
            max_expansions=retro_star_max_expansions,
            scscore_weights=scscore_weights if scscore_weights is not None else DEFAULT_SCSCORE_WEIGHTS,
        )
        return rs_planner.search(
            target_smiles,
            target_inchikey,
            max_steps=max_steps,
            route_backend=route_backend,
            offer_lookup=offer_lookup,
            lookup_backend=lookup_backend,
            template_backend=template_backend,
            offer_cache=cache,
            target_mol=target_mol,
            currency=currency,
            desalt_fallback=desalt_fallback,
            default_step_yield=default_step_yield,
            stoich_excess=stoich_excess,
            target_scale_g=target_scale_g,
            max_pack_multiple=max_pack_multiple,
            offer_mol_weights=offer_mol_weights,
        )

    routes: list[Route] = []

    if route_backend in {"auto", "aizynthfinder"} and aizynth_backend and aizynth_backend.is_ready():
        az_routes = aizynth_backend.search(target_smiles, target_inchikey, max_steps=max_steps)
        for r in az_routes:
            rebuilt = _build_route(
                target_smiles,
                target_inchikey,
                r.steps,
                r.leaves,
                r.method,
                cache,
                offer_lookup,
                score=r.score,
                currency=currency,
                desalt_fallback=desalt_fallback,
                default_step_yield=default_step_yield,
                stoich_excess=stoich_excess,
                target_scale_g=target_scale_g,
                max_pack_multiple=max_pack_multiple,
                offer_mol_weights=offer_mol_weights,
            )
            routes.append(rebuilt)
        has_solved_az = any(r.leaves and all(r.leaves_purchasable) for r in routes)
        if route_backend == "aizynthfinder" or (route_backend == "auto" and has_solved_az):
            routes.sort(key=_score_route_for_selection)
            return routes

    if route_backend in {"auto", "lookup"} and lookup_backend and lookup_backend.available:
        routes.extend(
            _expand_lookup_routes(
                target_smiles,
                target_inchikey,
                max_steps=max_steps,
                lookup_backend=lookup_backend,
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
        )
        if route_backend == "lookup":
            routes.sort(key=_score_route_for_selection)
            return routes

    if route_backend in {"auto", "template"}:
        tb = template_backend or TemplateFallbackBackend()
        effective_depth = min(max_steps, 2)
        method_label = (
            "template_fallback (heuristic RDKit SMARTS disconnections, max 2 steps; "
            "ignores stereochemistry and selectivity)"
        )

        step1_list = tb.one_step(target_smiles, mol=target_mol)
        solved_1step = False
        for st1 in step1_list:
            r1 = _build_route(
                target_smiles,
                target_inchikey,
                [st1],
                list(st1.reactants),
                method_label,
                cache,
                offer_lookup,
                currency=currency,
                desalt_fallback=desalt_fallback,
                default_step_yield=default_step_yield,
                stoich_excess=stoich_excess,
                target_scale_g=target_scale_g,
                max_pack_multiple=max_pack_multiple,
                offer_mol_weights=offer_mol_weights,
            )
            routes.append(r1)
            if r1.leaves and all(r1.leaves_purchasable):
                solved_1step = True

        if not solved_1step and effective_depth >= 2:
            for st1 in step1_list[:6]:
                missing_idx = [
                    idx
                    for idx, leaf in enumerate(st1.reactants)
                    if not _is_purchasable(
                        leaf, cache, offer_lookup, desalt_fallback=desalt_fallback
                    )
                ]
                if len(missing_idx) != 1:
                    continue
                m_idx = missing_idx[0]
                missing_leaf = st1.reactants[m_idx]
                missing_nm = normalize_smiles(missing_leaf)
                other_leaves = [r for i, r in enumerate(st1.reactants) if i != m_idx]
                for st2 in tb.one_step(
                    missing_leaf, mol=missing_nm.mol if missing_nm else None
                )[:6]:
                    all_leaves = other_leaves + list(st2.reactants)
                    r2 = _build_route(
                        target_smiles,
                        target_inchikey,
                        [st1, st2],
                        all_leaves,
                        method_label,
                        cache,
                        offer_lookup,
                        currency=currency,
                        desalt_fallback=desalt_fallback,
                        default_step_yield=default_step_yield,
                        stoich_excess=stoich_excess,
                        target_scale_g=target_scale_g,
                        max_pack_multiple=max_pack_multiple,
                        offer_mol_weights=offer_mol_weights,
                    )
                    routes.append(r2)
                    if r2.leaves and all(r2.leaves_purchasable):
                        break

    routes.sort(key=_score_route_for_selection)
    return routes


def _is_purchasable(
    smiles: str,
    offer_cache: dict[str, list[Offer]],
    offer_lookup: OfferLookupFn,
    *,
    desalt_fallback: bool = True,
) -> bool:
    nm = normalize_smiles(smiles)
    if nm is None:
        return False
    offs = _resolve_leaf_offers(nm, offer_cache, offer_lookup, desalt_fallback=desalt_fallback)
    return len(offs) > 0
