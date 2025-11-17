from __future__ import annotations

import math

import numpy as np
import pytest

from src.offline.segment_integrator import SolverConfig
from src.offline.stiff_ode import integrate_local_system


def _make_solver() -> SolverConfig:
    return SolverConfig(method="BDF", rtol=1e-6, atol=1e-9, max_step=0.5, seed=None)


def test_integrate_local_system_returns_initial_state_for_zero_span() -> None:
    solver = _make_solver()
    y0 = np.array([1.0, 2.0])

    result = integrate_local_system(lambda t, y: -y, y0, 0.5, 0.5, solver)

    assert np.allclose(result, y0)
    assert result is not y0


def test_integrate_local_system_matches_linear_system_solution() -> None:
    solver = _make_solver()
    fast_rate = 75.0
    slow_rate = 0.1
    span = 1.25
    y0 = np.array([2.0, 4.0])

    def rhs(_: float, y: np.ndarray) -> np.ndarray:
        return np.array(
            [
                -fast_rate * y[0],
                -slow_rate * y[1],
            ]
        )

    result = integrate_local_system(rhs, y0, 0.0, span, solver, max_internal_step_days=0.1)

    expected = np.array(
        [
            y0[0] * math.exp(-fast_rate * span),
            y0[1] * math.exp(-slow_rate * span),
        ]
    )
    assert result == pytest.approx(expected, rel=1e-6)
