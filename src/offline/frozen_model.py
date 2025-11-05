"""Backward-compatible entry points for the frozen model runtime."""

from __future__ import annotations

from .entities import CompiledExpression, DoseEntry, RuleEntry, ScenarioResult
from .simulation import EVENT_LOG_FIELDS, SEMANTICS_VERSION, simulate_frozen_model
from .snapshot import FrozenModel, _parse_trigger, load_frozen_model, sha256_file, snapshot_digest

__all__ = [
    "EVENT_LOG_FIELDS",
    "SEMANTICS_VERSION",
    "CompiledExpression",
    "DoseEntry",
    "RuleEntry",
    "ScenarioResult",
    "FrozenModel",
    "load_frozen_model",
    "simulate_frozen_model",
    "sha256_file",
    "snapshot_digest",
    "_parse_trigger",
]
