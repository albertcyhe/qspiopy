r"""Utilities for ingesting QSP-IO JSON parameter catalogues.

The MATLAB implementation of QSP-IO stores its physiological constants in
JSON files where each entry may either contain a literal value or provide an
expression referencing previously defined parameters.  The expressions follow
the original toolbox convention of using ``p(i)`` placeholders to refer to the
``i``\ th element of the accompanying ``derived_from`` list.

For the purposes of reproducing the published tutorial experiments we need a
small, self-contained loader that understands these semantics and exposes the
result as a regular Python mapping.  The helper is intentionally lightweight –
only a subset of mathematical functions are supported and expressions are
evaluated eagerly once all dependencies are available – but this is sufficient
for the parameter files bundled with the public repository as well as the
TNBC-specific catalogue referenced in the paper.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping


_SAFE_MATH_NAMESPACE: Dict[str, float] = {
    name: getattr(math, name)
    for name in dir(math)
    if not name.startswith("_")
}
_SAFE_MATH_NAMESPACE.update({
    "pi": math.pi,
    "e": math.e,
})


@dataclass(frozen=True)
class Parameter:
    """Container describing a resolved parameter entry."""

    name: str
    value: float
    units: str
    description: str
    source: str


class ParameterSet(Mapping[str, float]):
    """Mapping-like wrapper around a resolved parameter catalogue."""

    def __init__(self, parameters: Mapping[str, Parameter]):
        self._parameters = dict(parameters)

    def __getitem__(self, key: str) -> float:  # type: ignore[override]
        return self._parameters[key].value

    def __iter__(self) -> Iterator[str]:  # type: ignore[override]
        return iter(self._parameters)

    def __len__(self) -> int:  # type: ignore[override]
        return len(self._parameters)

    def metadata(self, name: str) -> Parameter:
        """Return the full :class:`Parameter` record for *name*.

        This is handy when downstream tooling needs access to the units or the
        textual provenance captured in the original JSON file.
        """

        return self._parameters[name]


def load_parameter_set(path: Path | str) -> ParameterSet:
    """Load a JSON parameter catalogue and resolve derived expressions.

    Parameters
    ----------
    path:
        Path to the JSON file.  The function accepts either a :class:`Path`
        instance or a raw string for convenience.
    """

    resolved: MutableMapping[str, Parameter] = {}
    with Path(path).open("r", encoding="utf8") as handle:
        raw_entries = json.load(handle)

    pending: List[MutableMapping[str, object]] = list(raw_entries)
    # Resolve parameters iteratively until all dependencies are satisfied.
    while pending:
        progress = False
        next_round: List[MutableMapping[str, object]] = []
        for entry in pending:
            name = str(entry["name"])
            if name in resolved:
                progress = True
                continue

            if entry.get("value") is not None:
                value = float(entry["value"])
                resolved[name] = Parameter(
                    name=name,
                    value=value,
                    units=str(entry.get("units", "")),
                    description=str(entry.get("description", "")),
                    source=str(entry.get("source", "")),
                )
                progress = True
                continue

            derived_from = entry.get("derived_from") or []
            expression = entry.get("expression")
            if not derived_from or not expression:
                raise ValueError(
                    f"Parameter {name} is missing both an explicit value and "
                    "a derivation recipe"
                )

            if all(dep in resolved for dep in derived_from):
                args = [resolved[dep].value for dep in derived_from]
                value = _evaluate_expression(str(expression), args)
                resolved[name] = Parameter(
                    name=name,
                    value=float(value),
                    units=str(entry.get("units", "")),
                    description=str(entry.get("description", "")),
                    source=str(entry.get("source", "")),
                )
                progress = True
            else:
                next_round.append(entry)

        if not progress:
            missing = {
                str(entry["name"]): [dep for dep in entry.get("derived_from", []) if dep not in resolved]
                for entry in next_round
            }
            raise ValueError(
                "Unable to resolve parameter dependencies: " + ", ".join(
                    f"{name} -> {deps}" for name, deps in missing.items()
                )
            )
        pending = next_round

    return ParameterSet(resolved)


def _evaluate_expression(expression: str, args: Iterable[float]) -> float:
    """Evaluate the MATLAB-style ``p(i)`` expression safely.

    The JSON catalogues use a compact syntax where ``p(1)`` resolves to the
    first dependency, ``p(2)`` to the second, and so on.  The evaluator exposes a
    helper callable ``p`` inside a restricted mathematical namespace and uses
    Python's :func:`eval` to compute the result.  All builtin access is removed
    to prevent arbitrary code execution.
    """

    values = list(float(v) for v in args)

    def _positional(index: float) -> float:
        idx = int(index) - 1
        if idx < 0 or idx >= len(values):
            raise IndexError(
                f"Expression requested p({index}) but only {len(values)} arguments available"
            )
        return values[idx]

    safe_locals: Dict[str, object] = dict(_SAFE_MATH_NAMESPACE)
    safe_locals["p"] = _positional
    safe_locals["abs"] = abs

    translated = expression.replace("^", "**")

    return float(
        eval(  # noqa: S307 - expression language is constrained via safe_locals
            translated,
            {"__builtins__": {}},
            safe_locals,
        )
    )


def combine_parameter_sets(paths: Iterable[Path | str]) -> ParameterSet:
    """Merge multiple parameter catalogues into a single :class:`ParameterSet`.

    Later entries take precedence when the same parameter name is defined
    multiple times.  This mirrors the MATLAB tooling where project-specific
    catalogues can override global defaults.
    """

    merged: Dict[str, Parameter] = {}
    for path in paths:
        subset = load_parameter_set(path)
        for name in subset:
            merged[name] = subset.metadata(name)
    return ParameterSet(merged)


def estimate_initial_tumour_cells(params: Mapping[str, float]) -> float:
    """Return the estimated number of cancer cells for the initial diameter."""

    diameter_cm = params.get("initial_tumour_diameter")
    if diameter_cm is None:
        raise KeyError("initial_tumour_diameter parameter missing")
    rho_cell = params.get("rho_cell")
    if rho_cell is None:
        raise KeyError("rho_cell parameter missing")

    volume_cm3 = (math.pi / 6.0) * float(diameter_cm) ** 3
    return float(volume_cm3 * float(rho_cell))


def tumour_volume_liters(cancer_cells: float, params: Mapping[str, float], *,
                         t_cells: float = 0.0, dead_cells: float = 0.0) -> float:
    """Compute tumour volume in litres using the article's compartment formula."""

    vol_cancer = params.get("vol_cell")
    vol_tcell = params.get("vol_Tcell")
    v_tmin = params.get("V_Tmin")
    if vol_cancer is None or vol_tcell is None or v_tmin is None:
        raise KeyError("vol_cell, vol_Tcell and V_Tmin parameters are required")

    # ``vol_cell`` and ``vol_Tcell`` are expressed in micrometre^3 / cell, while
    # ``V_Tmin`` is stored in microlitres.  Convert everything to litres so the
    # resulting volume is consistent with downstream reporting utilities.
    vol_cancer_l = float(vol_cancer) * 1e-15
    vol_tcell_l = float(vol_tcell) * 1e-15
    v_tmin_l = float(v_tmin) * 1e-6

    return (
        vol_cancer_l * (cancer_cells + dead_cells)
        + vol_tcell_l * t_cells
        + v_tmin_l
    )


def tumour_diameter_from_volume(volume_l: float) -> float:
    """Convert a spherical tumour volume (litres) to diameter in centimetres."""

    volume_cm3 = volume_l * 1e3
    radius_cm = ((3.0 * volume_cm3) / (4.0 * math.pi)) ** (1.0 / 3.0)
    return 2.0 * radius_cm

