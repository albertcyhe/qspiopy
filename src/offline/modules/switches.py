"""Lightweight module switches applied on top of frozen snapshots."""

from __future__ import annotations

from typing import Dict, Iterable

from ..snapshot import FrozenModel


def apply_parameter_overrides(model: FrozenModel, overrides: Dict[str, float] | None) -> FrozenModel:
    """Return the model after applying parameter overrides."""
    if not overrides:
        return model
    model.parameters = dict(model.parameters)
    model.parameters.update(overrides)
    return model


def disable_repeated_assignments(model: FrozenModel, targets: Iterable[str] | None) -> FrozenModel:
    """Remove repeated assignment rules targeting the provided symbols."""
    if not targets:
        return model
    block = {target.strip() for target in targets if target}
    if not block:
        return model
    model.repeated_assignment_order = [
        rule for rule in model.repeated_assignment_order if rule.target not in block
    ]
    return model
