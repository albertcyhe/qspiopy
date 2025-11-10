"""Target-volume initial condition generator for frozen snapshots."""

from __future__ import annotations

import fnmatch
import math
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .aliases import inject_output_aliases
from .segment_integrator import SolverConfig
from .snapshot import FrozenModel
import logging
import time

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
    max_wall_seconds: Optional[float] = None
    reset_policy: str = "cancer_only"  # {"all_zero","cancer_only","custom"}
    preserve_patterns: Tuple[str, ...] = ()
    solver: SolverConfig = field(
        default_factory=lambda: SolverConfig(method="BDF", rtol=1e-6, atol=1e-9, max_step=np.inf, seed=None)
    )


def generate_initial_conditions(model: FrozenModel, *, opts: ICOptions) -> Tuple[np.ndarray, Dict[str, float]]:
    """Integrate an untreated tumour until it reaches the target diameter."""

    seed_identifier = _resolve_species_identifier(model, opts.seed_species_id)
    state0 = model.initial_state().astype(float)
    model.apply_initial_assignments_to_state(state0)

    def _should_preserve(identifier: str, species_name: Optional[str]) -> bool:
        policy = (opts.reset_policy or "").lower()
        token = (species_name or identifier or "").lower()
        if policy == "all_zero":
            return False
        if policy == "cancer_only":
            return not token.startswith("c")
        if policy == "custom":
            patterns = opts.preserve_patterns or ()
            return any(fnmatch.fnmatch(token, pattern.lower()) for pattern in patterns)
        return False

    for ident, idx in model.dynamic_indices.items():
        entry = model.species_lookup.get(ident)
        if _should_preserve(ident, entry.name if entry else None):
            continue
        state0[idx] = 0.0

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

    start_wall = time.monotonic()
    last_logged = 0.0
    current_state = state0
    current_time = 0.0
    prev_ctx = _reconcile_context(model, current_state.copy())
    inject_output_aliases(prev_ctx)
    prev_diam = _diameter_cm_from_volume_l(float(prev_ctx.get("V_T", 0.0)))
    if prev_diam >= opts.target_diameter_cm:
        return current_state.copy(), dict(prev_ctx)

    while current_time < float(opts.max_days):
        if opts.max_wall_seconds and (time.monotonic() - start_wall) > opts.max_wall_seconds:
            raise RuntimeError(
                f"Failed to reach target diameter {opts.target_diameter_cm} cm within wall-time {opts.max_wall_seconds}s"
            )
        span_end = min(float(opts.max_days), current_time + 5.0)
        sol = solve_ivp(
            fun=model.rhs,
            t_span=(current_time, span_end),
            y0=current_state,
            method=opts.solver.method,
            rtol=opts.solver.rtol,
            atol=opts.solver.atol,
            max_step=opts.solver.max_step,
            events=[diameter_event],
            dense_output=True,
        )
        if not sol.success:
            raise RuntimeError(f"IC integration failed at t={current_time}: {sol.message}")
        prev_state = current_state
        prev_time = current_time
        current_state = np.asarray(sol.y[:, -1], dtype=float)
        current_time = span_end
        ctx = _reconcile_context(model, current_state.copy())
        inject_output_aliases(ctx)
        diam = _diameter_cm_from_volume_l(float(ctx.get("V_T", 0.0)))
        if diam - last_logged >= 0.05 or current_time >= float(opts.max_days):
            logging.getLogger(__name__).info(
                "ic_progress t=%.1f d=%.3f/%.3f cm",
                current_time,
                diam,
                opts.target_diameter_cm,
            )
            last_logged = diam
        if sol.t_events and sol.t_events[0].size:
            y_star = np.asarray(sol.y_events[0][-1], dtype=float)
            ctx_star = _reconcile_context(model, y_star.copy())
            inject_output_aliases(ctx_star)
            return y_star, ctx_star
        if prev_diam < opts.target_diameter_cm <= diam:
            if sol.sol is not None:
                t_lo, t_hi = prev_time, current_time
                for _ in range(40):
                    t_mid = 0.5 * (t_lo + t_hi)
                    y_mid = np.asarray(sol.sol(t_mid), dtype=float)
                    ctx_mid = _reconcile_context(model, y_mid.copy())
                    inject_output_aliases(ctx_mid)
                    diam_mid = _diameter_cm_from_volume_l(float(ctx_mid.get("V_T", 0.0)))
                    if diam_mid >= opts.target_diameter_cm:
                        return y_mid, ctx_mid
                    t_lo = t_mid
            return current_state.copy(), dict(ctx)
        prev_diam = diam

    final_ctx = _reconcile_context(model, current_state.copy())
    inject_output_aliases(final_ctx)
    raise RuntimeError(
        f"Failed to reach target diameter {opts.target_diameter_cm} cm within {opts.max_days} days "
        f"(last diameter = {_diameter_cm_from_volume_l(float(final_ctx.get('V_T', 0.0))):.3f} cm)"
    )
