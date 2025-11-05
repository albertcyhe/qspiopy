"""Core dataclasses shared across the offline runtime."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import sympy as sp


@dataclass(frozen=True)
class ScenarioResult:
    """Container describing the output of a frozen snapshot simulation."""

    time_days: np.ndarray
    cancer_cells: np.ndarray
    dead_cells: np.ndarray
    t_cells: np.ndarray
    tumour_volume_l: np.ndarray
    tumour_diameter_cm: np.ndarray
    pd1_occupancy: np.ndarray
    tcell_density_per_ul: np.ndarray
    extras: Dict[str, np.ndarray] = field(default_factory=dict)

    def to_frame(self) -> pd.DataFrame:
        data = {
            "time_days": self.time_days,
            "cancer_cells": self.cancer_cells,
            "dead_cells": self.dead_cells,
            "t_cells": self.t_cells,
            "tumour_volume_l": self.tumour_volume_l,
            "tumour_diameter_cm": self.tumour_diameter_cm,
            "pd1_occupancy": self.pd1_occupancy,
            "tcell_density_per_ul": self.tcell_density_per_ul,
        }
        for column, values in self.extras.items():
            data[column] = values
        return pd.DataFrame(data)


@dataclass(frozen=True)
class SpeciesEntry:
    identifier: str
    name: str
    compartment: str
    initial_amount: float
    constant: bool
    boundary_condition: bool
    nonnegative: bool
    initial_units: str
    units: str
    interpreted_dimension: str


@dataclass(frozen=True)
class CompiledExpression:
    tokens: Tuple[str, ...]
    func: object
    sympy_expr: Optional[sp.Expr] = None

    def evaluate(self, context: Dict[str, float]) -> float:
        result = self.evaluate_raw(context)
        return float(result)

    def evaluate_raw(self, context: Dict[str, float]):
        if not self.tokens:
            return self.func()
        values = [context[token] for token in self.tokens]
        return self.func(*values)


@dataclass(frozen=True)
class RuleEntry:
    index: int
    name: str
    rule_type: str
    target: str
    expression: str
    compiled: Optional[CompiledExpression]
    dependencies: Tuple[str, ...]
    safe_target: Optional[str]
    algebraic_solver: Optional[CompiledExpression]


@dataclass(frozen=True)
class ReactionEntry:
    name: str
    expression: str
    compiled: CompiledExpression


@dataclass(frozen=True)
class EventAssignment:
    target: str
    expression: str
    compiled: Optional[CompiledExpression]


@dataclass(frozen=True)
class EventEntry:
    index: int
    name: str
    trigger_expression: str
    trigger_compiled: CompiledExpression
    trigger_boolean_compiled: Optional[CompiledExpression]
    direction: float
    delay_expression: str
    delay_compiled: Optional[CompiledExpression]
    delay_type: str
    assignments: Tuple[EventAssignment, ...]

    @property
    def assignments_text(self) -> str:
        parts: List[str] = []
        for assignment in self.assignments:
            expr = assignment.expression or ""
            parts.append(f"{assignment.target}={expr}")
        return "; ".join(parts)


@dataclass(frozen=True)
class OdeEntry:
    identifier: str
    expression: str
    compiled: CompiledExpression


@dataclass(frozen=True)
class DoseEntry:
    index: int
    name: str
    dose_type: str
    target: str
    amount: float
    amount_units: str
    start_time: float
    interval: float
    repeat_count: int
    rate: Optional[float] = None
    rate_units: str = ""
    duration: Optional[float] = None


@dataclass(frozen=True, order=True)
class ScheduledDose:
    time: float
    priority: int
    dose: DoseEntry = field(compare=False)
    amount: float = field(compare=False)
