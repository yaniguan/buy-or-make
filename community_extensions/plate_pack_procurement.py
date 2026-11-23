"""Shared-Synthon Combinatorial Cart Optimizer (Idea #06).

Inspired by Enamine REAL / AutoDNA shared building-block plate synthesis workflows.
Selects a budgeted cart of shared in-stock synthons (`USD` catalog pack prices;
strictly excluding `null`-priced synthons) via budgeted marginal-coverage lookahead
and local swap refinement so purchasing `k` shared building blocks unlocks the
maximum number of 1-step target analogues under a strict `$B` cart budget.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SynthonPack:
    """In-stock reagent/building-block pack with strict catalog pack price."""

    synthon_id: str
    role: str  # e.g. "amine_core" or "acid_diversifier"
    pack_price_usd: float | None
    in_stock: bool = True

    def is_purchasable_priced(self) -> bool:
        return (
            self.in_stock
            and self.pack_price_usd is not None
            and math.isfinite(self.pack_price_usd)
            and self.pack_price_usd > 0.0
        )


@dataclass(frozen=True)
class TargetAnalogueSpec:
    """1-step target analogue unlocked when all `required_synthons` are in the cart."""

    target_id: str
    required_synthons: tuple[str, ...]
    priority_weight: float = 1.0


@dataclass(frozen=True)
class CartOptimizationResult:
    """Selected building-block cart and unlocked 1-step target analogue summary."""

    purchased_synthon_ids: tuple[str, ...]
    unlocked_target_ids: tuple[str, ...]
    total_spend_usd: float
    cart_budget_usd: float
    targets_unlocked_per_100_usd: float
    synthon_reuse_factor: float


def _count_unlocked(
    cart_ids: set[str],
    targets: Sequence[TargetAnalogueSpec],
) -> list[TargetAnalogueSpec]:
    return [t for t in targets if all(sid in cart_ids for sid in t.required_synthons)]


def optimize_independent_per_target_cart(
    synthons: Sequence[SynthonPack],
    targets: Sequence[TargetAnalogueSpec],
    cart_budget_usd: float = 500.0,
) -> CartOptimizationResult:
    """Baseline buying missing synthons target-by-target without cross-library lookahead."""
    if cart_budget_usd <= 0.0:
        raise ValueError("cart_budget_usd must be positive")
    priced_map = {s.synthon_id: s for s in synthons if s.is_purchasable_priced()}

    cart: set[str] = set()
    spend = 0.0

    for target in targets:
        if any(sid not in priced_map for sid in target.required_synthons):
            continue
        missing = [sid for sid in target.required_synthons if sid not in cart]
        inc_cost = sum(float(priced_map[sid].pack_price_usd or 0.0) for sid in missing)
        if spend + inc_cost <= cart_budget_usd + 1e-9:
            cart.update(missing)
            spend += inc_cost

    unlocked = _count_unlocked(cart, targets)
    per_100 = len(unlocked) / (cart_budget_usd / 100.0)
    reuse = sum(len(t.required_synthons) for t in unlocked) / len(cart) if cart else 0.0
    return CartOptimizationResult(
        purchased_synthon_ids=tuple(sorted(cart)),
        unlocked_target_ids=tuple(t.target_id for t in unlocked),
        total_spend_usd=round(spend, 2),
        cart_budget_usd=cart_budget_usd,
        targets_unlocked_per_100_usd=round(per_100, 4),
        synthon_reuse_factor=round(reuse, 4),
    )


def optimize_shared_synthon_cart(
    synthons: Sequence[SynthonPack],
    targets: Sequence[TargetAnalogueSpec],
    cart_budget_usd: float = 500.0,
) -> CartOptimizationResult:
    """Budgeted combinatorial maximum-coverage cart optimizer over shared synthons.

    Combines 1-synthon and 2-synthon marginal unlock-per-dollar greedy expansion
    (handling zero-marginal cold start of bimolecular 1-step reactions) with
    2-step degree potential tie-breaking and budget-feasible swap refinement.
    """
    if cart_budget_usd <= 0.0:
        raise ValueError("cart_budget_usd must be positive")

    priced_map = {s.synthon_id: s for s in synthons if s.is_purchasable_priced()}
    feasible_targets = [
        t for t in targets if all(sid in priced_map for sid in t.required_synthons)
    ]

    # Degree potential: how many feasible library targets reference each synthon
    degree: dict[str, float] = {sid: 0.0 for sid in priced_map}
    for t in feasible_targets:
        for sid in t.required_synthons:
            degree[sid] += t.priority_weight

    cart: set[str] = set()
    spend = 0.0

    while True:
        remaining = cart_budget_usd - spend
        cur_unlocked_count = len(_count_unlocked(cart, feasible_targets))

        best_bundle: tuple[str, ...] | None = None
        best_density = 0.0
        best_potential = -1.0

        # Evaluate singletons for marginal immediate unlocks once cores are present
        for sid, pack in priced_map.items():
            if sid in cart:
                continue
            cost = float(pack.pack_price_usd or 0.0)
            if cost > remaining + 1e-9:
                continue
            trial = cart | {sid}
            gain = len(_count_unlocked(trial, feasible_targets)) - cur_unlocked_count
            if gain > 0:
                density = gain / cost
                pot = degree[sid] / cost
                if (density, pot) > (best_density, best_potential):
                    best_density = density
                    best_potential = pot
                    best_bundle = (sid,)

        # Evaluate missing synthon bundles per target (handles 2-synthon cold start!)
        for t in feasible_targets:
            missing = tuple(
                sorted(sid for sid in t.required_synthons if sid not in cart)
            )
            if not missing:
                continue
            cost = sum(float(priced_map[sid].pack_price_usd or 0.0) for sid in missing)
            if cost > remaining + 1e-9:
                continue
            trial = cart | set(missing)
            gain = len(_count_unlocked(trial, feasible_targets)) - cur_unlocked_count
            # Include fractional combinatorial hub bonus on cold-start pairs
            hub_potential = sum(degree[sid] for sid in missing) / cost
            density = gain / cost
            if (density, hub_potential) > (best_density, best_potential):
                best_density = density
                best_potential = hub_potential
                best_bundle = missing

        if best_bundle is None:
            break

        bundle_cost = sum(
            float(priced_map[sid].pack_price_usd or 0.0) for sid in best_bundle
        )
        cart.update(best_bundle)
        spend += bundle_cost

    # Local 1-for-1 swap refinement under strict cart budget
    improved = True
    while improved:
        improved = False
        cur_count = len(_count_unlocked(cart, feasible_targets))
        for drop_id in sorted(cart):
            drop_cost = float(priced_map[drop_id].pack_price_usd or 0.0)
            for add_id, add_pack in priced_map.items():
                if add_id in cart:
                    continue
                add_cost = float(add_pack.pack_price_usd or 0.0)
                new_spend = spend - drop_cost + add_cost
                if new_spend > cart_budget_usd + 1e-9:
                    continue
                trial_cart = (cart - {drop_id}) | {add_id}
                new_count = len(_count_unlocked(trial_cart, feasible_targets))
                if new_count > cur_count or (
                    new_count == cur_count and new_spend < spend - 1e-6
                ):
                    cart = trial_cart
                    spend = new_spend
                    improved = True
                    break
            if improved:
                break

    unlocked = _count_unlocked(cart, feasible_targets)
    per_100 = len(unlocked) / (cart_budget_usd / 100.0)
    reuse = sum(len(t.required_synthons) for t in unlocked) / len(cart) if cart else 0.0
    return CartOptimizationResult(
        purchased_synthon_ids=tuple(sorted(cart)),
        unlocked_target_ids=tuple(t.target_id for t in unlocked),
        total_spend_usd=round(spend, 2),
        cart_budget_usd=cart_budget_usd,
        targets_unlocked_per_100_usd=round(per_100, 4),
        synthon_reuse_factor=round(reuse, 4),
    )


def run_benchmark(cart_budget_usd: float = 500.0, seed: int = 42) -> dict[str, object]:
    """Compare `independent_per_target_buy` vs `shared_synthon_cart` under a `$500` budget."""
    rng = np.random.default_rng(seed)

    # 10 amine cores ($42-$58/pack) + 10 carboxylic acid diversifiers ($38-$55/pack)
    # + 12 bespoke single-use specialty synthons ($75-$95/pack) + 3 unpriced synthons (`None`)
    synthons: list[SynthonPack] = []
    for i in range(10):
        synthons.append(
            SynthonPack(
                f"AMINE_CORE_{i}",
                "amine_core",
                round(float(rng.uniform(42.0, 56.0)), 2),
            )
        )
    for j in range(10):
        synthons.append(
            SynthonPack(
                f"ACID_DIV_{j}",
                "acid_diversifier",
                round(float(rng.uniform(38.0, 54.0)), 2),
            )
        )
    for k in range(12):
        synthons.append(
            SynthonPack(
                f"BESPOKE_BB_{k}",
                "bespoke_singleton",
                round(float(rng.uniform(76.0, 92.0)), 2),
            )
        )
    for u in range(3):
        synthons.append(
            SynthonPack(f"UNPRICED_BB_{u}", "unpriced", None, in_stock=True)
        )

    targets: list[TargetAnalogueSpec] = []
    # Incoming request queue starts with bespoke non-shared analogues (trap for naive per-target buyer!)
    for k in range(0, 12, 2):
        targets.append(
            TargetAnalogueSpec(
                target_id=f"T_bespoke_{k // 2}",
                required_synthons=(f"BESPOKE_BB_{k}", f"BESPOKE_BB_{k + 1}"),
            )
        )

    # Dense combinatorial array over shared Enamine-style cores x diversifiers (64 1-step analogues)
    for i in range(8):
        for j in range(8):
            if (i + j) % 5 == 4:
                continue
            targets.append(
                TargetAnalogueSpec(
                    target_id=f"T_grid_{i}_{j}",
                    required_synthons=(f"AMINE_CORE_{i}", f"ACID_DIV_{j}"),
                )
            )

    # Unpriced target analogues (must never be hallucinated as purchasable!)
    for u in range(3):
        targets.append(
            TargetAnalogueSpec(
                target_id=f"T_unpriced_{u}",
                required_synthons=(f"AMINE_CORE_{u}", f"UNPRICED_BB_{u}"),
            )
        )

    indep = optimize_independent_per_target_cart(
        synthons, targets, cart_budget_usd=cart_budget_usd
    )
    shared = optimize_shared_synthon_cart(
        synthons, targets, cart_budget_usd=cart_budget_usd
    )

    return {
        "idea_id": "idea_06_plate_pack_procurement",
        "cart_budget_usd": cart_budget_usd,
        "total_candidate_targets": len(targets),
        "independent_per_target_buy": {
            "target_analogues_unlocked": len(indep.unlocked_target_ids),
            "target_analogues_per_100_usd": indep.targets_unlocked_per_100_usd,
            "total_spend_usd": indep.total_spend_usd,
            "synthons_purchased": len(indep.purchased_synthon_ids),
            "synthon_reuse_factor": indep.synthon_reuse_factor,
        },
        "shared_synthon_cart": {
            "target_analogues_unlocked": len(shared.unlocked_target_ids),
            "target_analogues_per_100_usd": shared.targets_unlocked_per_100_usd,
            "total_spend_usd": shared.total_spend_usd,
            "synthons_purchased": len(shared.purchased_synthon_ids),
            "synthon_reuse_factor": shared.synthon_reuse_factor,
            "budget_respected": shared.total_spend_usd <= cart_budget_usd,
        },
    }
