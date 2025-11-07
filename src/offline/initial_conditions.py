"""Target-volume initial condition generator for frozen snapshots."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .segment_integrator import SolverConfig
from .snapshot import FrozenModel

CM3_PER_LITRE = 1000.0  # 1 L = 1000 cm^3


def _diameter_cm_from_volume_l(volume_l: float) -> float:
    if volume_l <= 0.0:
        return 0.0
    volume_cm3 = volume_l * CM3_PER_LITRE
    radius = ((3.0 * volume_cm3) / (4.0 * math.pi)) ** (1.0 / 3.0)
    return 2.0 * radius


def _reconcile_context(model: FrozenModel, vec: np.ndarray) -> Dict[str, float]:
    ctx = model.build_context_from_state(vec.copy())
    model.evaluate_repeated_assignments(ctx)
    model.apply_algebraic_rules(ctx, vec, mutate=False)
    model.sync_state_from_context(ctx, vec)
    return ctx


def _resolve_species_identifier(model: FrozenModel, token: str) -> str:
    if token in model.dynamic_indices:
        return token
    entry = model.species_lookup.get(token)
    if entry and entry.identifier in model.dynamic_indices:
        return entry.identifier
    entry = model.species_name_lookup.get(token)
    if entry and entry.identifier in model.dynamic_indices:
        return entry.identifier
    raise ValueError(f"seed species '{token}' not found in dynamic indices")


@dataclass(frozen=True)
class ICOptions:
    """Options controlling the target-volume IC generation."""

    target_diameter_cm: float
    seed_species_id: str = "C1"
    seed_cells: float = 1.0
    max_days: float = 4000.0
    solver: SolverConfig = field(
        default_factory=lambda: SolverConfig(method="BDF", rtol=1e-6, atol=1e-9, max_step=np.inf, seed=None)
    )


def generate_initial_conditions(model: FrozenModel, *, opts: ICOptions) -> Tuple[np.ndarray, Dict[str, float]]:
    """Integrate an untreated tumour until it reaches the target diameter."""

    seed_identifier = _resolve_species_identifier(model, opts.seed_species_id)
    state0 = model.initial_state().astype(float)
    model.apply_initial_assignments_to_state(state0)
    state0[:] = 0.0
    seed_idx = model.dynamic_indices[seed_identifier]
    state0[seed_idx] = float(opts.seed_cells)

    def diameter_event(time: float, values: np.ndarray) -> float:
        vec = np.asarray(values, dtype=float)
        ctx = model.build_context_from_state(vec.copy())
        model.evaluate_repeated_assignments(ctx)
        model.apply_algebraic_rules(ctx, vec, mutate=False)
        volume_l = float(ctx.get("V_T", 0.0))
        return _diameter_cm_from_volume_l(volume_l) - opts.target_diameter_cm

    diameter_event.direction = 1.0
    diameter_event.terminal = True

    sol = solve_ivp(
        fun=model.rhs,
        t_span=(0.0, float(opts.max_days)),
        y0=state0,
        method=opts.solver.method,
        rtol=opts.solver.rtol,
        atol=opts.solver.atol,
        max_step=opts.solver.max_step,
        events=[diameter_event],
        dense_output=False,
    )

    if not sol.t_events or not sol.t_events[0].size:
        last_ctx = _reconcile_context(model, sol.y[:, -1].copy())
        raise RuntimeError(
            f"Failed to reach target diameter {opts.target_diameter_cm} cm within {opts.max_days} days "
            f"(last diameter = {_diameter_cm_from_volume_l(float(last_ctx.get('V_T', 0.0))):.3f} cm)"
        )

    y_star = np.asarray(sol.y_events[0][-1], dtype=float)
    ctx_star = _reconcile_context(model, y_star.copy())
    return y_star, ctx_star
