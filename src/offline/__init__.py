"""Stable public API for the frozen SimBiology runtime."""

from .entities import DoseEntry, ScenarioResult
from .simulation import simulate_frozen_model
from .snapshot import load_frozen_model

__all__ = (
    "load_frozen_model",
    "simulate_frozen_model",
    "ScenarioResult",
    "DoseEntry",
)
