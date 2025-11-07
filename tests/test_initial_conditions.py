from __future__ import annotations

import numpy as np
import pytest

from src.offline.initial_conditions import ICOptions, _diameter_cm_from_volume_l, generate_initial_conditions
from src.offline.segment_integrator import SolverConfig


class DummySpecies:
    def __init__(self, identifier: str):
        self.identifier = identifier
        self.name = identifier
        self.constant = False
        self.boundary_condition = False


class DummyModel:
    def __init__(self, growth_rate: float, cell_volume_l: float):
        self.growth_rate = growth_rate
        self.cell_volume_l = cell_volume_l
        self.dynamic_indices = {"C1": 0}
        self.species_lookup = {"C1": DummySpecies("C1")}
        self.species_name_lookup = {"C1": DummySpecies("C1")}
        self.events: list = []
        self.doses: list = []
        self.config: dict = {}

    def initial_state(self) -> np.ndarray:
        return np.array([0.0], dtype=float)

    def apply_initial_assignments_to_state(self, state: np.ndarray) -> None:  # pragma: no cover - no-op
        return

    def build_context_from_state(self, state: np.ndarray):
        value = float(state[0])
        return {"C1": value, "V_T": value * self.cell_volume_l}

    def evaluate_repeated_assignments(self, context):  # pragma: no cover - no-op
        return

    def apply_algebraic_rules(self, context, state, mutate=False):  # pragma: no cover - no-op
        return

    def sync_state_from_context(self, context, state):
        state[0] = context.get("C1", state[0])

    def rhs(self, t, y):
        return np.array([self.growth_rate * y[0]], dtype=float)


@pytest.mark.slow
def test_generate_initial_conditions_hits_target():
    model = DummyModel(growth_rate=0.05, cell_volume_l=1e-12)
    opts = ICOptions(
        target_diameter_cm=0.02,
        seed_cells=1.0,
        solver=SolverConfig(method="BDF", rtol=1e-8, atol=1e-12, max_step=0.5, seed=None),
    )
    state, context = generate_initial_conditions(model, opts=opts)
    diameter = _diameter_cm_from_volume_l(float(context.get("V_T", 0.0)))
    assert abs(diameter - opts.target_diameter_cm) < 1e-3
    assert state[0] > 0.0
