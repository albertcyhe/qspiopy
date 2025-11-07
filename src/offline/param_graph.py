"""Parameter dependency graph for derived SimBiology parameters."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Mapping, Optional, Sequence

from .errors import SnapshotError


def _evaluate_expression(expr: str, values: Sequence[float]) -> float:
    substituted = expr or ""
    for idx, val in enumerate(values, start=1):
        substituted = substituted.replace(f"p({idx})", f"({val})")
    try:
        safe_globals = {"__builtins__": None, "pi": math.pi, "e": math.e}
        return float(eval(substituted, safe_globals, {}))
    except Exception as exc:  # pragma: no cover - expressions mirror snapshot exports
        raise SnapshotError(f"Failed to evaluate expression '{expr}': {exc}") from exc


@dataclass
class ParameterSpec:
    name: str
    unit: str = ""
    expression: str = ""
    dependencies: Sequence[str] = field(default_factory=list)
    value: Optional[float] = None


class ParameterGraph:
    """Manages base and derived parameters with dependency resolution."""

    def __init__(self, converter: Callable[[float, str], float]):
        self._convert = converter
        self._values: Dict[str, float] = {}
        self._metadata: Dict[str, Dict[str, object]] = {}
        self._specs: Dict[str, ParameterSpec] = {}

    def add_base(self, name: str, value: float, unit: str) -> None:
        canonical = self._convert(value, unit)
        self._values[name] = canonical
        self._metadata[name] = {
            "type": "base",
            "unit": unit,
            "raw_value": value,
        }

    def add_spec(
        self,
        name: str,
        unit: str = "",
        *,
        expression: str = "",
        dependencies: Optional[Sequence[str]] = None,
        value: Optional[float] = None,
    ) -> None:
        self._specs[name] = ParameterSpec(
            name=name,
            unit=unit,
            expression=expression or "",
            dependencies=list(dependencies or []),
            value=value,
        )

    def evaluate(self) -> Mapping[str, float]:
        remaining = dict(self._specs)
        while remaining:
            progress = False
            for name, spec in list(remaining.items()):
                deps = spec.dependencies
                if any(dep not in self._values for dep in deps):
                    continue
                if spec.expression:
                    dep_values = [self._values[dep] for dep in deps]
                    raw_value = _evaluate_expression(spec.expression, dep_values)
                elif spec.value is not None:
                    raw_value = float(spec.value)
                else:
                    raise SnapshotError(f"Derived parameter '{name}' lacks expression/value")
                canonical = self._convert(raw_value, spec.unit)
                self._values[name] = canonical
                self._metadata[name] = {
                    "type": "derived",
                    "unit": spec.unit,
                    "raw_value": raw_value,
                    "expression": spec.expression,
                    "dependencies": list(spec.dependencies),
                }
                remaining.pop(name)
                progress = True
            if not progress:
                unresolved = ", ".join(sorted(remaining.keys()))
                raise SnapshotError(f"Unable to resolve derived parameters: {unresolved}")
        return dict(self._values)

    @property
    def metadata(self) -> Mapping[str, Dict[str, object]]:
        return self._metadata
