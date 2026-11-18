"""Retrosynthetic route backends: reaction lookup (USPTO/ORD), template fallback, optional AiZynthFinder."""

from __future__ import annotations

import json
import os
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .canonical import normalize_smiles
from .cost import EST_COST_METHOD, estimate_route_cost
from .schema import Offer, Route, RouteStep

DEFAULT_REACTION_DB = Path("data/reactions.sqlite")

# Curated open SMARTS retrosynthesis disconnections (product >> reactants).
# Limits: heuristic structural disconnections only; ignores stereochemistry, protecting
# groups, regioselectivity, chemoselectivity, and reaction conditions.
RETRO_TEMPLATES: list[tuple[str, str]] = [
    (
        "amide_coupling_secondary",
        "[#6:1][C:2](=[O:3])[NH1:4][#6:5]>>[#6:1][C:2](=[O:3])O.[NH2:4][#6:5]",
    ),
    (
        "amide_coupling_tertiary",
        "[#6:1][C:2](=[O:3])[N:4]([#6:5])[#6:6]>>[#6:1][C:2](=[O:3])O.[NH1:4]([#6:5])[#6:6]",
    ),
    (
        "esterification",
        "[#6:1][C:2](=[O:3])[O:4][#6:5]>>[#6:1][C:2](=[O:3])O.[OH1:4][#6:5]",
    ),
    (
        "sulfonamide_secondary",
        "[#6:1][S:2](=[O:3])(=[O:4])[NH1:5][#6:6]>>[#6:1][S:2](=[O:3])(=[O:4])Cl.[NH2:5][#6:6]",
    ),
    (
        "sulfonamide_tertiary",
        "[#6:1][S:2](=[O:3])(=[O:4])[N:5]([#6:6])[#6:7]>>[#6:1][S:2](=[O:3])(=[O:4])Cl.[NH1:5]([#6:6])[#6:7]",
    ),
    (
        "urea_formation",
        "[#6:1][NH1:2][C:3](=[O:4])[NH1:5][#6:6]>>[#6:1][N:2]=[C:3]=[O:4].[NH2:5][#6:6]",
    ),
    (
        "carbamate_formation",
        "[#6:1][O:2][C:3](=[O:4])[NH1:5][#6:6]>>[#6:1][OH1:2].[#6:6][N:5]=[C:3]=[O:4]",
    ),
    (
        "suzuki_biaryl",
        "[c:1]-[c:2]>>[c:1]Br.[c:2]B(O)O",
    ),
    (
        "suzuki_aryl_alkenyl",
        "[c:1]-[CX3:2]=[CX3:3]>>[c:1]Br.[CX3:2](B(O)O)=[CX3:3]",
    ),
    (
        "williamson_aryl_ether",
        "[c:1][O:2][CX4:3]>>[c:1][OH1:2].[CX4:3]Br",
    ),
    (
        "williamson_alkyl_ether",
        "[CX4:1][O:2][CX4:3]>>[CX4:1][OH1:2].[CX4:3]Br",
    ),
    (
        "reductive_amination_secondary",
        "[#6:1][CH2:2][NH1:3][#6:4]>>[#6:1][CH1:2]=O.[NH2:3][#6:4]",
    ),
    (
        "reductive_amination_tertiary",
        "[#6:1][CH2:2][N:3]([#6:4])[#6:5]>>[#6:1][CH1:2]=O.[NH1:3]([#6:4])[#6:5]",
    ),
    (
        "buchwald_hartwig_cn",
        "[c:1][NH1:2][c:3]>>[c:1]Br.[NH2:2][c:3]",
    ),
]


@dataclass(frozen=True)
class _CompiledTemplate:
    name: str
    smarts: str
    rxn: rdChemReactions.ChemicalReaction


_COMPILED_TEMPLATES: Optional[list[_CompiledTemplate]] = None


def _get_compiled_templates() -> list[_CompiledTemplate]:
    global _COMPILED_TEMPLATES
    if _COMPILED_TEMPLATES is None:
        _COMPILED_TEMPLATES = [
            _CompiledTemplate(name=n, smarts=s, rxn=rdChemReactions.ReactionFromSmarts(s))
            for n, s in RETRO_TEMPLATES
        ]
    return _COMPILED_TEMPLATES


class ReactionLookupBackend:
    """Exact product-InChIKey lookup in locally indexed USPTO (CC0) / ORD (CC-BY-SA-4.0) reactions."""

    def __init__(self, db_path: Path | str = DEFAULT_REACTION_DB) -> None:
        self.db_path = Path(db_path)
        self._conn: Optional[sqlite3.Connection] = None
        if self.db_path.exists():
            self._conn = sqlite3.connect(str(self.db_path))

    @property
    def available(self) -> bool:
        return self._conn is not None

    def one_step(self, inchikey: str, canonical_smiles: str, limit: int = 10) -> list[RouteStep]:
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
        if self._conn is not None:
            self._conn.close()
            self._conn = None


class TemplateFallbackBackend:
    """Built-in 1- and 2-step heuristic SMARTS retrosynthetic disconnection engine."""

    def __init__(self, max_branch_per_node: int = 12) -> None:
        self.max_branch_per_node = max_branch_per_node
        self.templates = _get_compiled_templates()

    def one_step(self, canonical_smiles: str) -> list[RouteStep]:
        mol = Chem.MolFromSmiles(canonical_smiles)
        if mol is None:
            return []
        parent_heavy = mol.GetNumHeavyAtoms()
        candidates: list[tuple[int, RouteStep]] = []
        seen_sets: set[tuple[str, ...]] = set()

        for tpl in self.templates:
            try:
                outcomes = tpl.rxn.RunReactants((mol,))
            except Exception:
                continue
            for prod_tuple in outcomes:
                reactants: list[str] = []
                total_heavy = 0
                ok = True
                for frag in prod_tuple:
                    try:
                        Chem.SanitizeMol(frag)
                        smi = Chem.MolToSmiles(frag, isomericSmiles=True, canonical=True)
                    except Exception:
                        ok = False
                        break
                    nm = normalize_smiles(smi)
                    if nm is None or nm.canonical_smiles == canonical_smiles:
                        ok = False
                        break
                    m_frag = Chem.MolFromSmiles(nm.canonical_smiles)
                    if m_frag is None or m_frag.GetNumHeavyAtoms() < 2:
                        ok = False
                        break
                    total_heavy += m_frag.GetNumHeavyAtoms()
                    reactants.append(nm.canonical_smiles)

                if not ok or not reactants:
                    continue
                # Prevent degenerate expansions that blow up size without progress
                if any(Chem.MolFromSmiles(r).GetNumHeavyAtoms() >= parent_heavy + 3 for r in reactants):
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
                # Prefer balanced disconnections (smaller max fragment heavy atoms)
                max_frag = max(Chem.MolFromSmiles(r).GetNumHeavyAtoms() for r in reactants)
                candidates.append((max_frag, step))

        candidates.sort(key=lambda item: (item[0], len(item[1].reactants)))
        return [step for _, step in candidates[: self.max_branch_per_node]]


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
            self._finder.expansion_policy.select_all()
            self._finder.filter_policy.select_all()
            self._finder.stock.select_all()

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
            routes.append(
                Route(
                    target_inchikey=target_inchikey,
                    n_steps=len(steps),
                    steps=steps,
                    leaves=leaves,
                    leaves_purchasable=leaves_ok,
                    method="aizynthfinder (MCTS retrosynthesis)",
                    score=route_dict.get("score"),
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
    if node.get("type") == "mol":
        children = node.get("children", [])
        if not children:
            leaves.append(node["smiles"])
            leaves_ok.append(bool(node.get("in_stock", False)))
            return
        rxn_node = children[0]
        reactant_nodes = rxn_node.get("children", [])
        reactant_smiles = [c["smiles"] for c in reactant_nodes if c.get("type") == "mol"]
        steps.append(
            RouteStep(
                product=node["smiles"],
                reactants=reactant_smiles,
                reaction_smiles=f"{'.'.join(reactant_smiles)}>>{node['smiles']}",
                template_or_reference=str(rxn_node.get("metadata", {}).get("template_code", "aizynthfinder")),
                source="aizynthfinder",
            )
        )
        for child in reactant_nodes:
            _walk_aizynth_tree(child, steps, leaves, leaves_ok)


OfferLookupFn = Callable[[str, str], list[Offer]]


def _score_route_for_selection(route: Route) -> tuple:
    all_ok = 0 if (route.leaves and all(route.leaves_purchasable)) else 1
    cost_tier = 2
    cost_val = 0.0
    if route.est_cost is not None and route.est_cost.value is not None:
        cost_tier = 0 if route.est_cost_complete else 1
        cost_val = route.est_cost.value
    return (all_ok, route.n_steps, cost_tier, cost_val, len(route.leaves))


def _build_route(
    target_smiles: str,
    target_inchikey: str,
    steps: list[RouteStep],
    leaves: list[str],
    method: str,
    offer_cache: dict[str, list[Offer]],
    offer_lookup: OfferLookupFn,
) -> Route:
    leaves_purchasable: list[bool] = []
    leaf_offers_by_ik: dict[str, list[Offer]] = {}
    for leaf_smi in leaves:
        nm = normalize_smiles(leaf_smi)
        if nm is None:
            leaves_purchasable.append(False)
            continue
        if nm.inchikey not in offer_cache:
            offer_cache[nm.inchikey] = offer_lookup(nm.inchikey, nm.canonical_smiles)
        offs = offer_cache[nm.inchikey]
        leaf_offers_by_ik[nm.inchikey] = offs
        leaves_purchasable.append(len(offs) > 0)

    est_cost, est_complete = estimate_route_cost(target_smiles, leaves, leaf_offers_by_ik)
    return Route(
        target_inchikey=target_inchikey,
        n_steps=len(steps),
        steps=steps,
        leaves=leaves,
        leaves_purchasable=leaves_purchasable,
        method=method,
        score=None,
        est_cost=est_cost,
        est_cost_method=EST_COST_METHOD,
        est_cost_complete=est_complete,
    )


def find_routes(
    target_smiles: str,
    target_inchikey: str,
    *,
    max_steps: int = 3,
    route_backend: str = "auto",
    offer_lookup: OfferLookupFn,
    lookup_backend: Optional[ReactionLookupBackend] = None,
    template_backend: Optional[TemplateFallbackBackend] = None,
    aizynth_backend: Optional[AiZynthFinderBackend] = None,
    offer_cache: Optional[dict[str, list[Offer]]] = None,
) -> list[Route]:
    """Search for retrosynthetic routes up to max_steps using the selected backend."""
    if max_steps <= 0:
        return []

    cache = offer_cache if offer_cache is not None else {}
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
            )
            routes.append(rebuilt)
        if routes or route_backend == "aizynthfinder":
            routes.sort(key=_score_route_for_selection)
            return routes

    if route_backend in {"auto", "lookup"} and lookup_backend and lookup_backend.available:
        l_steps = lookup_backend.one_step(target_inchikey, target_smiles)
        for st in l_steps:
            routes.append(
                _build_route(
                    target_smiles,
                    target_inchikey,
                    [st],
                    list(st.reactants),
                    f"reaction_lookup ({st.source})",
                    cache,
                    offer_lookup,
                )
            )
        if route_backend == "lookup":
            routes.sort(key=_score_route_for_selection)
            return routes

    if route_backend in {"auto", "template"}:
        tb = template_backend or TemplateFallbackBackend()
        # Cap heuristic SMARTS expansion at 2 steps to stay fast and avoid combinatorial noise
        effective_depth = min(max_steps, 2)
        method_label = (
            "template_fallback (heuristic RDKit SMARTS disconnections, max 2 steps; "
            "ignores stereochemistry and selectivity)"
        )

        step1_list = tb.one_step(target_smiles)
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
            )
            routes.append(r1)
            if r1.leaves and all(r1.leaves_purchasable):
                solved_1step = True

        # If not solved in 1 step and depth >= 2, try disconnecting non-purchasable leaves
        if not solved_1step and effective_depth >= 2:
            for st1 in step1_list[:6]:
                # Check which reactants of st1 are not purchasable
                missing_idx = [
                    idx
                    for idx, leaf in enumerate(st1.reactants)
                    if not _is_purchasable(leaf, cache, offer_lookup)
                ]
                # Only expand if exactly 1 reactant is missing (keeps 2-step search crisp and clean)
                if len(missing_idx) != 1:
                    continue
                m_idx = missing_idx[0]
                missing_leaf = st1.reactants[m_idx]
                other_leaves = [r for i, r in enumerate(st1.reactants) if i != m_idx]
                for st2 in tb.one_step(missing_leaf)[:6]:
                    all_leaves = other_leaves + list(st2.reactants)
                    r2 = _build_route(
                        target_smiles,
                        target_inchikey,
                        [st1, st2],
                        all_leaves,
                        method_label,
                        cache,
                        offer_lookup,
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
) -> bool:
    nm = normalize_smiles(smiles)
    if nm is None:
        return False
    if nm.inchikey not in offer_cache:
        offer_cache[nm.inchikey] = offer_lookup(nm.inchikey, nm.canonical_smiles)
    return len(offer_cache[nm.inchikey]) > 0
