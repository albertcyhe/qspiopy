"""Shared stiff ODE helpers reused by white-box modules."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Callable, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .errors import NumericsError
if TYPE_CHECKING:  # pragma: no cover
    from .segment_integrator import SolverConfig

StateVector = np.ndarray
RhsFn = Callable[[float, StateVector], StateVector]
EventFns = Optional[Sequence[Callable[[float, np.ndarray], float]]]


def _looks_like_step_failure(message: str) -> bool:
    text = (message or "").lower()
    return ("step size" in text) or ("strictly increasing" in text)


def solve_stiff_ivp(
    rhs: RhsFn,
    span: Tuple[float, float],
    y0: StateVector,
    solver: "SolverConfig",
    *,
    first_step: Optional[float] = None,
    max_step: Optional[float] = None,
    dense_output: bool = False,
    events: EventFns = None,
    jac_sparsity: Optional[np.ndarray] = None,
    vectorized: bool = False,
    allow_shrink: bool = True,
    max_attempts: int = 8,
):
    """Wrapper around solve_ivp with consistent retry semantics."""

    t0 = float(span[0])
    t1 = float(span[1])
    if t1 == t0:
        dummy = solve_ivp(
            rhs,
            (t0, t1),
            np.asarray(y0, dtype=float),
            method=solver.method,
            rtol=solver.rtol,
            atol=solver.atol,
            max_step=1.0,
        )
        return dummy

    state0 = np.asarray(y0, dtype=float)
    total_span = abs(t1 - t0)
    attempt_first = None if first_step is None or first_step <= 0.0 else float(first_step)

    attempt_max = max_step
    if attempt_max is None or attempt_max <= 0.0 or not math.isfinite(attempt_max):
        cap = float(solver.max_step or 0.0)
        if cap > 0.0 and math.isfinite(cap):
            attempt_max = cap
        else:
            attempt_max = total_span
    min_cap = max(total_span * 1e-6, 1e-12)
    result = None
    for _ in range(max_attempts):
        result = solve_ivp(
            rhs,
            (t0, t1),
            state0,
            method=solver.method,
            rtol=solver.rtol,
            atol=solver.atol,
            max_step=attempt_max,
            first_step=attempt_first,
            dense_output=dense_output,
            events=events,
            jac_sparsity=jac_sparsity,
            vectorized=vectorized,
        )
        if result.success or not allow_shrink:
            return result
        if not _looks_like_step_failure(result.message or ""):
            return result
        attempt_max = max(attempt_max * 0.5, min_cap)
        if attempt_first is not None:
            attempt_first = min(attempt_first, attempt_max)
    return result if result is not None else solve_ivp(
        rhs,
        (t0, t1),
        state0,
        method=solver.method,
        rtol=solver.rtol,
        atol=solver.atol,
        max_step=attempt_max,
        first_step=attempt_first,
        dense_output=dense_output,
        events=events,
        jac_sparsity=jac_sparsity,
        vectorized=vectorized,
    )


def integrate_local_system(
    rhs: RhsFn,
    y0: StateVector,
    t0: float,
    t1: float,
    solver: "SolverConfig",
    *,
    max_internal_step_days: Optional[float] = 1e-4,
) -> StateVector:
    """Integrate a small stiff system y' = rhs(t, y) reusing SolverConfig."""
    start = float(t0)
    stop = float(t1)
    span = stop - start
    state0 = np.array(y0, dtype=float, copy=True)
    if span <= 0.0 or state0.size == 0:
        return state0

    max_step = max_internal_step_days
    if max_step is not None:
        if max_step <= 0.0:
            max_step = None
        else:
            solver_cap = float(solver.max_step or 0.0)
            if solver_cap > 0.0 and math.isfinite(solver_cap):
                max_step = min(max_step, solver_cap)

    sol = solve_stiff_ivp(
        rhs,
        (start, stop),
        state0,
        solver,
        max_step=max_step,
        dense_output=False,
        events=None,
        allow_shrink=True,
    )
    if not sol.success or not sol.y.size:
        raise NumericsError(f"Stiff ODE integration failed: {sol.message}")
    return np.asarray(sol.y[:, -1], dtype=float)


__all__ = ["integrate_local_system", "solve_stiff_ivp"]
