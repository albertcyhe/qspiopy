"""Public exports for the frozen SimBiology runtime."""

from .frozen_model import ScenarioResult, simulate_frozen_model

__all__ = [
    "ScenarioResult",
    "simulate_frozen_model",
]
