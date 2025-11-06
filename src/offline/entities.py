"""Core dataclasses shared across the offline runtime."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import sympy as sp


SEMANTICS_VERSION = "1.0"
BASE_HEADER = (
    "time_days",
    "cancer_cells",
    "dead_cells",
    "t_cells",
    "tumour_volume_l",
    "tumour_diameter_cm",
    "pd1_occupancy",
    "tcell_density_per_ul",
)


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
    semantics_version: str = SEMANTICS_VERSION
    header: Tuple[str, ...] = field(default_factory=tuple)
    provenance: Dict[str, str] = field(default_factory=dict)

    def column_order(self) -> Tuple[str, ...]:
        if self.header:
            return self.header
        extras = tuple(self.extras.keys())
        return BASE_HEADER + extras

    def to_frame(self, order: str = "contract") -> pd.DataFrame:
        columns = self.column_order() if order == "contract" else tuple(
            list(BASE_HEADER) + list(self.extras.keys())
        )
        data_map = {
            "time_days": self.time_days,
            "cancer_cells": self.cancer_cells,
            "dead_cells": self.dead_cells,
            "t_cells": self.t_cells,
            "tumour_volume_l": self.tumour_volume_l,
            "tumour_diameter_cm": self.tumour_diameter_cm,
            "pd1_occupancy": self.pd1_occupancy,
            "tcell_density_per_ul": self.tcell_density_per_ul,
        }
        data_map.update(self.extras)
        frame = pd.DataFrame({column: data_map[column] for column in columns})
        frame.attrs["semantics_version"] = self.semantics_version
        if self.provenance:
            frame.attrs["provenance"] = self.provenance
        return frame

    def save_csv(
        self,
        path: Path,
        *,
        order: str = "contract",
        include_header_manifest: bool = True,
        **to_csv_kwargs,
    ) -> None:
        frame = self.to_frame(order=order)
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False, **to_csv_kwargs)
        if include_header_manifest:
            manifest_path = path.with_suffix(path.suffix + ".header.txt")
            manifest_path.write_text("\n".join(self.column_order()), encoding="utf8")


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
    amount_mg: Optional[float] = None


@dataclass(frozen=True, order=True)
class ScheduledDose:
    time: float
    priority: int
    dose: DoseEntry = field(compare=False)
    amount: float = field(compare=False)
    amount_mg: Optional[float] = field(default=None, compare=False)
