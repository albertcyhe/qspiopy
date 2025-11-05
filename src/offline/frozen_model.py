"""Load and simulate frozen SimBiology snapshots without MATLAB."""

from __future__ import annotations

import hashlib
import heapq
import json
import logging
import math
from collections import deque
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import sympy as sp

SAFE_GLOBALS = {
    "__builtins__": None,
    "exp": math.exp,
    "log": math.log,
    "log10": math.log10,
    "sqrt": math.sqrt,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "sinh": math.sinh,
    "cosh": math.cosh,
    "tanh": math.tanh,
    "abs": abs,
    "min": min,
    "max": max,
    "pow": pow,
    "pi": math.pi,
}

SEMANTICS_VERSION = "0.9.0"

EVENT_LOG_FIELDS = (
    "event_index",
    "time_fire",
    "time_trigger",
    "delay",
    "type",
    "assignments",
)

_TIME_TOL = 1e-9
_T0_REL_TOL = 1e-9
_T0_ABS_TOL = 1e-12


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


@dataclass(frozen=True, order=True)
class ScheduledEvent:
    time: float
    priority: int
    entry: EventEntry = field(compare=False)
    trigger_time: float = field(compare=False)
    delay: float = field(compare=False)
    assignments: str = field(compare=False)


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
    dose: DoseEntry
    amount: float


@dataclass
class FrozenModel:
    name: str
    source_dir: Path
    species: List[SpeciesEntry]
    species_lookup: Dict[str, SpeciesEntry]
    species_name_lookup: Dict[str, SpeciesEntry]
    parameters: Dict[str, float]
    compartments: Dict[str, float]
    config: Dict[str, object]
    rules: List[RuleEntry]
    rate_rules: List[RuleEntry]
    algebraic_rules: List[RuleEntry]
    reactions: List[ReactionEntry]
    odes: List[OdeEntry]
    events: List[EventEntry]
    doses: List[DoseEntry]
    variants: List[Dict[str, object]]
    stoichiometry: List[List[Tuple[int, float]]]
    constants: Dict[str, float]
    dynamic_indices: Dict[str, int]
    repeated_assignment_order: List[RuleEntry]
    initial_assignment_rules: List[RuleEntry]

    def initial_state(self) -> np.ndarray:
        state = np.zeros(len(self.dynamic_indices), dtype=float)
        for species in self.species:
            if species.identifier in self.dynamic_indices:
                state[self.dynamic_indices[species.identifier]] = species.initial_amount
        return state

    def _base_context(self) -> Dict[str, float]:
        ctx = {}
        ctx.update(self.parameters)
        ctx.update(self.compartments)
        ctx.update(self.constants)
        return ctx

    def build_context_from_state(self, state: np.ndarray) -> Dict[str, float]:
        context = self._base_context()
        for identifier, idx in self.dynamic_indices.items():
            value = float(state[idx])
            context[identifier] = value
            entry = self.species_lookup.get(identifier)
            if entry is not None:
                context[entry.name] = value
        return context

    def apply_initial_assignments(self, context: Dict[str, float]) -> None:
        for rule in self.initial_assignment_rules:
            context[rule.target] = rule.compiled.evaluate(context)

    def evaluate_repeated_assignments(self, context: Dict[str, float]) -> None:
        for rule in self.repeated_assignment_order:
            context[rule.target] = rule.compiled.evaluate(context)

    def apply_algebraic_rules(
        self,
        context: Dict[str, float],
        state: Optional[np.ndarray],
        *,
        mutate: bool = False,
        atol: float = 1e-9,
    ) -> None:
        if not self.algebraic_rules:
            return
        for rule in self.algebraic_rules:
            if rule.algebraic_solver is not None:
                raw_value = rule.algebraic_solver.evaluate(context)
                value = float(raw_value)
                if mutate:
                    self._apply_target_value(rule.target, value, context, state)
                else:
                    self._set_context_value(rule.target, value, context, state)
                continue
            if rule.compiled is None:
                continue
            residual = float(rule.compiled.evaluate(context))
            if abs(residual) > atol:
                raise RuntimeError(
                    f"Algebraic rule '{rule.name or rule.target}' residual {residual} exceeds tolerance {atol}"
                )

    def _set_context_value(
        self,
        target: str,
        value: float,
        context: Dict[str, float],
        state: Optional[np.ndarray] = None,
    ) -> None:
        if target in self.parameters or target in self.compartments:
            context[target] = value
            return
        if target in self.dynamic_indices:
            idx = self.dynamic_indices[target]
            if state is not None:
                state[idx] = value
            context[target] = value
            entry = self.species_lookup.get(target)
            if entry is not None:
                context[entry.name] = value
            return
        if target in self.species_name_lookup:
            entry = self.species_name_lookup[target]
            identifier = entry.identifier
            if identifier in self.dynamic_indices:
                idx = self.dynamic_indices[identifier]
                if state is not None:
                    state[idx] = value
            context[identifier] = value
            context[target] = value
            return
        context[target] = value

    def _apply_target_value(
        self, target: str, value: float, context: Dict[str, float], state: Optional[np.ndarray] = None
    ) -> None:
        if target in self.parameters:
            self.parameters[target] = value
            context[target] = value
            return
        if target in self.compartments:
            self.compartments[target] = value
            context[target] = value
            return
        if target in self.dynamic_indices:
            idx = self.dynamic_indices[target]
            if state is not None:
                state[idx] = value
            context[target] = value
            entry = self.species_lookup.get(target)
            if entry is not None:
                context[entry.name] = value
            return
        if target in self.species_name_lookup:
            entry = self.species_name_lookup[target]
            identifier = entry.identifier
            if identifier in self.dynamic_indices:
                idx = self.dynamic_indices[identifier]
                if state is not None:
                    state[idx] = value
            context[identifier] = value
            context[target] = value
            return
        context[target] = value

    def apply_initial_assignments_to_state(self, state: np.ndarray) -> None:
        if not self.initial_assignment_rules:
            return
        context = self.build_context_from_state(state.copy())
        for rule in self.initial_assignment_rules:
            value = rule.compiled.evaluate(context)
            self._apply_target_value(rule.target, value, context, state)
        self.enforce_species_constraints(context, state)

    def sync_state_from_context(self, context: Dict[str, float], state: np.ndarray) -> None:
        for identifier, idx in self.dynamic_indices.items():
            if identifier in context:
                state[idx] = context[identifier]
                name_entry = self.species_lookup.get(identifier)
                if name_entry is not None:
                    context[name_entry.name] = context[identifier]
            else:
                name_entry = self.species_lookup.get(identifier)
                if name_entry is not None and name_entry.name in context:
                    state[idx] = context[name_entry.name]
        self.enforce_species_constraints(context, state)

    def enforce_species_constraints(self, context: Dict[str, float], state: np.ndarray) -> None:
        for identifier, idx in self.dynamic_indices.items():
            entry = self.species_lookup.get(identifier)
            if entry is None:
                continue
            if entry.nonnegative and state[idx] < 0.0:
                state[idx] = 0.0
                context[identifier] = 0.0
                context[entry.name] = 0.0

    def evaluate_reactions(self, context: Dict[str, float]) -> Dict[str, float]:
        values: Dict[str, float] = {}
        for reaction in self.reactions:
            try:
                value = reaction.compiled.evaluate(context)
            except Exception as exc:  # pragma: no cover
                raise RuntimeError(f"Failed to evaluate reaction {reaction.name}: {reaction.expression}") from exc
            context[reaction.name] = value
            values[reaction.name] = value
        return values

    def evaluate_ode_rhs(self, context: Dict[str, float]) -> np.ndarray:
        derivative = np.zeros(len(self.dynamic_indices), dtype=float)
        for ode in self.odes:
            target_idx = self.dynamic_indices.get(ode.identifier)
            if target_idx is None:
                continue
            entry = self.species_lookup.get(ode.identifier)
            if entry is not None and entry.boundary_condition:
                derivative[target_idx] = 0.0
                continue
            derivative[target_idx] = ode.compiled.evaluate(context)
        for rule in self.rate_rules:
            if rule.compiled is None:
                continue
            idx = self.dynamic_indices.get(rule.target)
            if idx is None:
                continue
            derivative[idx] = rule.compiled.evaluate(context)
        return derivative

    def apply_dose(self, dose: DoseEntry, amount: float, context: Dict[str, float], state: np.ndarray) -> None:
        target = dose.target
        entry = self.species_lookup.get(target)
        if entry is None and target in self.species_name_lookup:
            entry = self.species_name_lookup[target]
        if entry is None:
            new_value = context.get(target, 0.0) + amount
            self._apply_target_value(target, new_value, context, state)
            return
        delta = amount
        units_lower = entry.units.lower()
        dimension_lower = entry.interpreted_dimension.lower()
        if (
            "molar" in units_lower
            or "concentration" in dimension_lower
            or "amount/vol" in dimension_lower
            or "/l" in units_lower
        ):
            compartment_name = entry.compartment
            compartment_volume = context.get(compartment_name, self.compartments.get(compartment_name))
            if compartment_volume in (None, 0.0):
                raise RuntimeError(f"Missing compartment volume for dose target {target}")
            delta = amount / compartment_volume
        current_value = context.get(entry.identifier, 0.0)
        new_value = current_value + delta
        self._apply_target_value(entry.identifier, new_value, context, state)

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        state_view = np.array(y, dtype=float, copy=True)
        context = self.build_context_from_state(state_view)
        self.evaluate_repeated_assignments(context)
        self.apply_algebraic_rules(context, state_view, mutate=False)
        self.evaluate_reactions(context)
        derivative = self.evaluate_ode_rhs(context)
        return derivative


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def _load_flux_expressions(path: Path) -> Dict[str, str]:
    if not path.is_file():
        return {}
    with path.open("r", encoding="utf8") as handle:
        lines = [line.strip() for line in handle]
    fluxes: Dict[str, str] = {}
    in_flux_section = False
    for line in lines:
        if not line:
            if in_flux_section:
                break
            continue
        if line.startswith("Fluxes:"):
            in_flux_section = True
            continue
        if in_flux_section:
            if "=" in line:
                name, expr = line.split("=", 1)
                fluxes[name.strip()] = expr.strip()
            else:
                break
    return fluxes


def _load_configset(path: Path) -> Dict[str, object]:
    if not path.is_file():
        return {}
    with path.open("r", encoding="utf8") as handle:
        return json.load(handle)


def _sanitize_tokens(tokens: Iterable[str]) -> List[str]:
    seen: Dict[str, str] = {}
    ordered: List[str] = []
    for token in tokens:
        if token not in seen and token:
            seen[token] = token
            ordered.append(token)
    ordered.sort(key=len, reverse=True)
    return ordered


def _load_ode_expressions(path: Path) -> List[Tuple[str, str]]:
    if not path.is_file():
        return []
    lines = [line.rstrip() for line in path.read_text(encoding="utf8").splitlines()]
    odes: List[Tuple[str, str]] = []
    in_section = False
    for line in lines:
        if not line.strip():
            if in_section:
                break
            continue
        if line.startswith("ODEs:"):
            in_section = True
            continue
        if line.startswith("Fluxes:"):
            break
        if in_section and line.startswith("d(") and ")/dt" in line:
            left, rhs = line.split("=", 1)
            identifier = left[2 : left.index(")/dt")].strip()
            odes.append((identifier, rhs.strip()))
    return odes


def _compile_expression(
    expression: str,
    tokens: Sequence[str],
    safe_names: Dict[str, str],
    reverse_safe: Dict[str, str],
) -> CompiledExpression:
    raw = expression.replace("^", "**")
    temp = raw
    placeholders: Dict[str, str] = {}
    for idx, token in enumerate(tokens):
        placeholder = f"__PLACEHOLDER_{idx}__"
        if token not in safe_names:
            safe_name = f"SYM_{len(safe_names)}"
            safe_names[token] = safe_name
            reverse_safe[safe_name] = token
        temp = temp.replace(token, placeholder)
        placeholders[placeholder] = safe_names[token]
    for placeholder, name in placeholders.items():
        temp = temp.replace(placeholder, name)
    sym_expr = sp.sympify(temp)
    free_symbols = sorted(sym_expr.free_symbols, key=lambda s: s.name)
    tokens_order = tuple(reverse_safe[str(sym)] for sym in free_symbols)
    func = sp.lambdify(free_symbols, sym_expr, modules=["math"])
    return CompiledExpression(tokens=tokens_order, func=func, sympy_expr=sym_expr)


def _compile_constant_expression(expression: str) -> CompiledExpression:
    expr = sp.sympify(expression.replace("^", "**"))
    func = sp.lambdify([], expr, modules=["math"])
    return CompiledExpression(tokens=(), func=func, sympy_expr=expr)


def _load_species(path: Path) -> Tuple[List[SpeciesEntry], Dict[str, float], List[str]]:
    frame = _read_csv(path)
    species: List[SpeciesEntry] = []
    constants: Dict[str, float] = {}
    for row in frame.itertuples(index=False):
        constant = str(row.constant_amount).lower() == "true"
        boundary = str(row.boundary_condition).lower() == "true"
        nonnegative = str(row.nonnegative).lower() == "true"
        initial_units = ""
        if hasattr(row, "initial_units") and not pd.isna(row.initial_units):
            initial_units = str(row.initial_units)
        units = ""
        if hasattr(row, "units") and not pd.isna(row.units):
            units = str(row.units)
        interpreted = ""
        if hasattr(row, "interpreted_dimension") and not pd.isna(row.interpreted_dimension):
            interpreted = str(row.interpreted_dimension)
        entry = SpeciesEntry(
            identifier=str(row.identifier),
            name=str(row.name),
            compartment=str(row.compartment),
            initial_amount=float(row.initial_amount),
            constant=constant,
            boundary_condition=boundary,
            nonnegative=nonnegative,
            initial_units=initial_units,
            units=units,
            interpreted_dimension=interpreted,
        )
        species.append(entry)
        if constant:
            constants[entry.identifier] = entry.initial_amount
            constants[entry.name] = entry.initial_amount
    tokens = []
    for entry in species:
        tokens.append(entry.identifier)
        tokens.append(entry.name)
    return species, constants, tokens


def _load_parameters(path: Path) -> Tuple[Dict[str, float], List[str]]:
    frame = _read_csv(path)
    params: Dict[str, float] = {}
    for row in frame.itertuples(index=False):
        name = str(row.name)
        value = float(row.value)
        units = str(row.units) if hasattr(row, "units") else ""
        params[name] = _convert_parameter_value(value, units)
    return params, list(params.keys())


def _load_compartments(path: Path) -> Tuple[Dict[str, float], List[str]]:
    frame = _read_csv(path)
    comps: Dict[str, float] = {}
    for row in frame.itertuples(index=False):
        name = str(row.name)
        value = float(row.capacity)
        units = str(row.units) if hasattr(row, "units") else ""
        comps[name] = _convert_compartment_value(value, units)
    return comps, list(comps.keys())


def _convert_compartment_value(value: float, units: str) -> float:
    if not units:
        return value
    u = units.lower()
    if u == "microliter":
        return value * 1e-6
    if u == "millimeter^3":
        return value * 1e-6
    if u == "litre" or u == "liter":
        return value
    return value


def _convert_parameter_value(value: float, units: str) -> float:
    if not units:
        return value
    u = units.lower()
    if u == "dimensionless":
        return value
    if u == "microliter":
        return value * 1e-6
    if u == "millimeter^3":
        return value * 1e-6
    if u == "micrometer^3":
        return value * 1e-15
    if u == "micrometer^3/cell":
        return value * 1e-15
    if u == "litre" or u == "liter":
        return value
    if u == "cell/milliliter":
        return value * 1000.0
    if u == "micromolarity":
        return value * 1e-6
    if u == "nanomolarity":
        return value * 1e-9
    if u == "nanomole/cell/hour":
        return value * 1e-9 * 24.0
    if u == "1/minute":
        return value * 1440.0
    if u == "1/second":
        return value * 86400.0
    if u == "1/(centimeter^3*minute)":
        return value * 1e3 * 1440.0
    if u == "1/day/milliliter":
        return value * 1000.0
    if u == "1/(molarity*second)":
        return value * 86400.0
    if u == "1/(micromolarity*nanometer*second)":
        return value * (1e6) * (1e-9) * 86400.0
    return value


def _load_rules(
    path: Path,
    tokens: Sequence[str],
    safe_names: Dict[str, str],
    reverse_safe: Dict[str, str],
) -> List[RuleEntry]:
    frame = _read_csv(path)
    rules: List[RuleEntry] = []
    for row in frame.itertuples(index=False):
        rule_type = str(row.type)
        target = "" if pd.isna(row.target) else str(row.target)
        raw_expr = "" if pd.isna(row.expression) else str(row.expression)
        rule_type_lower = rule_type.strip().lower()
        compiled: Optional[CompiledExpression] = None
        if raw_expr:
            if rule_type_lower in {
                "repeatedassignment",
                "initialassignment",
                "rate",
                "raterule",
                "algebraic",
                "algebraicrule",
            }:
                compiled = _compile_expression(raw_expr, tokens, safe_names, reverse_safe)
        dep_tokens: Tuple[str, ...] = tuple()
        if compiled is not None:
            dep_tokens = tuple(token for token in compiled.tokens if token != target)
        safe_target = safe_names.get(target) if target else None
        algebraic_solver: Optional[CompiledExpression] = None
        if rule_type_lower in {"algebraic", "algebraicrule"} and compiled is not None and safe_target:
            try:
                target_symbol = sp.Symbol(safe_target)
                solutions = sp.solve(sp.Eq(compiled.sympy_expr, 0), target_symbol, dict=True)
            except Exception:  # pragma: no cover - best effort
                solutions = []
            if solutions:
                solution_expr = solutions[0][target_symbol]
                free_symbols = sorted(solution_expr.free_symbols, key=lambda s: s.name)
                solver_tokens = tuple(reverse_safe[str(sym)] for sym in free_symbols)
                solver_func = sp.lambdify(free_symbols, solution_expr, modules=["math"])
                algebraic_solver = CompiledExpression(
                    tokens=solver_tokens,
                    func=solver_func,
                    sympy_expr=solution_expr,
                )
        rules.append(
            RuleEntry(
                index=int(row.rule_index),
                name=str(row.name),
                rule_type=rule_type,
                target=target,
                expression=raw_expr,
                compiled=compiled,
                dependencies=dep_tokens,
                safe_target=safe_target,
                algebraic_solver=algebraic_solver,
            )
        )
    return rules


def _parse_trigger(trigger: str) -> Tuple[str, float]:
    operators = ["<=", ">=", "<", ">", "=="]
    for op in operators:
        if op in trigger:
            left, right = trigger.split(op, 1)
            left = left.strip()
            right = right.strip()
            if op in ("<", "<="):
                direction = -1.0
            elif op in (">", ">="):
                direction = 1.0
            else:
                direction = 0.0
            expr = f"({left})-({right})"
            return expr, direction
    return trigger, 0.0


def _parse_assignments(assignments: str) -> List[Tuple[str, str]]:
    if not assignments:
        return []
    pairs: List[Tuple[str, str]] = []
    for part in assignments.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            target, expr = part.split("=", 1)
            pairs.append((target.strip(), expr.strip()))
    return pairs


def _load_events(
    path: Path,
    tokens: Sequence[str],
    safe_names: Dict[str, str],
    reverse_safe: Dict[str, str],
) -> List[EventEntry]:
    if not path.is_file():
        return []
    frame = _read_csv(path)
    entries: List[EventEntry] = []
    for row in frame.itertuples(index=False):
        trigger_expr_raw = str(row.trigger)
        expr_str, direction = _parse_trigger(trigger_expr_raw)
        trigger_compiled = _compile_expression(expr_str, tokens, safe_names, reverse_safe)
        trigger_bool_compiled: Optional[CompiledExpression] = None
        if trigger_expr_raw:
            trigger_bool_compiled = _compile_expression(trigger_expr_raw, tokens, safe_names, reverse_safe)
        delay_expr = "" if pd.isna(row.delay) else str(row.delay)
        delay_compiled = _compile_expression(delay_expr, tokens, safe_names, reverse_safe) if delay_expr else None
        delay_type_raw = "" if not hasattr(row, "delay_type") or pd.isna(row.delay_type) else str(row.delay_type)
        delay_type = delay_type_raw.strip().lower()
        assignments_raw = _parse_assignments(str(row.assignments))
        assignments: List[EventAssignment] = []
        for target, expr in assignments_raw:
            if expr:
                compiled = _compile_expression(expr, tokens, safe_names, reverse_safe)
            else:
                compiled = None
            assignments.append(EventAssignment(target=target, expression=expr, compiled=compiled))
        entries.append(
            EventEntry(
                index=int(row.event_index),
                name=str(row.name),
                trigger_expression=trigger_expr_raw,
                trigger_compiled=trigger_compiled,
                trigger_boolean_compiled=trigger_bool_compiled,
                direction=direction,
                delay_expression=delay_expr,
                delay_compiled=delay_compiled,
                delay_type=delay_type,
                assignments=tuple(assignments),
            )
        )
    return entries


def _load_doses(path: Path) -> List[DoseEntry]:
    if not path.is_file():
        return []
    frame = _read_csv(path)
    doses: List[DoseEntry] = []
    for row in frame.itertuples(index=False):
        amount = 0.0 if pd.isna(row.amount) else float(row.amount)
        amount_units = "" if pd.isna(row.amount_units) else str(row.amount_units)
        rate = None
        if hasattr(row, "rate") and not pd.isna(row.rate):
            rate = float(row.rate)
        rate_units = ""
        if hasattr(row, "rate_units") and not pd.isna(row.rate_units):
            rate_units = str(row.rate_units)
        duration = None
        if hasattr(row, "duration") and not pd.isna(row.duration):
            duration = float(row.duration)
        doses.append(
            DoseEntry(
                index=int(row.dose_index),
                name=str(row.name),
                dose_type=str(row.type),
                target=str(row.target),
                amount=amount,
                amount_units=amount_units,
                start_time=float(row.start_time),
                interval=float(row.interval),
                repeat_count=int(row.repeat_count),
                rate=rate,
                rate_units=rate_units,
                duration=duration,
            )
        )
    return doses


def _load_variants(path: Path) -> List[Dict[str, object]]:
    if not path.is_file():
        return []
    frame = _read_csv(path)
    variants: List[Dict[str, object]] = []
    for row in frame.itertuples(index=False):
        variants.append(
            {
                "variant_index": int(row.variant_index),
                "name": str(row.name),
                "active": str(row.active).lower() == "true",
                "content": str(row.content),
            }
        )
    return variants


def _load_reactions(
    path: Path,
    tokens: Sequence[str],
    safe_names: Dict[str, str],
    reverse_safe: Dict[str, str],
    fluxes: Dict[str, str],
) -> List[ReactionEntry]:
    frame = _read_csv(path)
    reactions: List[ReactionEntry] = []
    for row in frame.itertuples(index=False):
        expr = fluxes.get(str(row.name), str(row.reaction_rate))
        if expr:
            compiled = _compile_expression(expr, tokens, safe_names, reverse_safe)
        else:
            compiled = _compile_constant_expression("0")
        reactions.append(ReactionEntry(name=str(row.name), expression=expr, compiled=compiled))
    return reactions


def _topological_order_repeated(rules: List[RuleEntry]) -> List[RuleEntry]:
    repeated = [rule for rule in rules if rule.rule_type.lower() == "repeatedassignment" and rule.compiled is not None]
    targets = {rule.target: rule for rule in repeated}
    in_degree: Dict[str, int] = {}
    adjacency: Dict[str, List[str]] = {target: [] for target in targets}
    for rule in repeated:
        deps = [dep for dep in rule.dependencies if dep in targets]
        in_degree[rule.target] = len(deps)
        for dep in deps:
            adjacency.setdefault(dep, []).append(rule.target)
    queue = deque([target for target, deg in in_degree.items() if deg == 0])
    order_targets: List[str] = []
    while queue:
        current = queue.popleft()
        order_targets.append(current)
        for neighbour in adjacency.get(current, []):
            in_degree[neighbour] -= 1
            if in_degree[neighbour] == 0:
                queue.append(neighbour)
    if len(order_targets) != len(repeated):
        remaining = [target for target, deg in in_degree.items() if deg > 0]
        raise ValueError(f"Cycle detected in repeated assignments: {remaining}")
    return [targets[target] for target in order_targets]


def _load_stoichiometry(path: Path, reaction_order: List[str], dynamic_indices: Dict[str, int]) -> List[List[Tuple[int, float]]]:
    frame = _read_csv(path)
    contributions: List[List[Tuple[int, float]]] = [[] for _ in reaction_order]
    for row in frame.itertuples(index=False):
        species = str(row.species)
        reaction = str(row.reaction)
        value = float(row.stoichiometry)
        if species not in dynamic_indices:
            continue
        if reaction not in reaction_order:
            continue
        reaction_idx = reaction_order.index(reaction)
        species_idx = dynamic_indices[species]
        contributions[reaction_idx].append((species_idx, value))
    return contributions


@lru_cache(maxsize=None)
def load_frozen_model(name: str, *, root: Path = Path("artifacts") / "matlab_frozen_model") -> FrozenModel:
    model_dir = root / name
    if not model_dir.is_dir():
        raise FileNotFoundError(f"Snapshot directory {model_dir} is missing. Run export_matlab_snapshot first.")

    species, constant_species, species_tokens = _load_species(model_dir / "species.csv")
    parameters, parameter_tokens = _load_parameters(model_dir / "parameters.csv")
    compartments, compartment_tokens = _load_compartments(model_dir / "compartments.csv")
    config = _load_configset(model_dir / "configset.json")

    all_tokens = _sanitize_tokens(
        list(species_tokens) + parameter_tokens + compartment_tokens + list(constant_species.keys())
    )

    tokens_sorted = _sanitize_tokens(all_tokens)
    safe_names = {token: f"SYM_{index}" for index, token in enumerate(tokens_sorted)}
    reverse_safe = {safe: token for token, safe in safe_names.items()}

    fluxes = _load_flux_expressions(model_dir / "equations.txt")

    rules = _load_rules(model_dir / "rules.csv", tokens_sorted, safe_names, reverse_safe)
    reactions = _load_reactions(model_dir / "reactions.csv", tokens_sorted, safe_names, reverse_safe, fluxes)
    reaction_tokens = [reaction.name for reaction in reactions]
    ode_defs = _load_ode_expressions(model_dir / "equations.txt")
    ode_tokens = _sanitize_tokens(tokens_sorted + reaction_tokens)
    odes: List[OdeEntry] = []
    for identifier, expr in ode_defs:
        compiled = _compile_expression(expr, ode_tokens, safe_names, reverse_safe)
        odes.append(OdeEntry(identifier=identifier, expression=expr, compiled=compiled))

    dynamic_species = [entry for entry in species if not entry.constant]
    dynamic_indices = {entry.identifier: idx for idx, entry in enumerate(dynamic_species)}
    stoichiometry = _load_stoichiometry(model_dir / "stoichiometry.csv", [r.name for r in reactions], dynamic_indices)

    constants = constant_species.copy()

    species_lookup = {entry.identifier: entry for entry in species}
    species_name_lookup = {entry.name: entry for entry in species}
    events = _load_events(model_dir / "events.csv", ode_tokens, safe_names, reverse_safe)
    doses = _load_doses(model_dir / "doses.csv")
    variants = _load_variants(model_dir / "variants.csv")
    repeated_order = _topological_order_repeated(rules)
    initial_rules = [rule for rule in rules if rule.rule_type.lower() == "initialassignment" and rule.compiled is not None]
    rate_rules = [rule for rule in rules if rule.rule_type.lower() in {"rate", "raterule"} and rule.compiled is not None]
    algebraic_rules = [
        rule for rule in rules if rule.rule_type.lower() in {"algebraic", "algebraicrule"} and rule.compiled is not None
    ]

    return FrozenModel(
        name=name,
        source_dir=model_dir,
        species=species,
        species_lookup=species_lookup,
        species_name_lookup=species_name_lookup,
        parameters=parameters,
        compartments=compartments,
        config=config,
        rules=rules,
        rate_rules=rate_rules,
        algebraic_rules=algebraic_rules,
        reactions=reactions,
        odes=odes,
        events=events,
        doses=doses,
        variants=variants,
        stoichiometry=stoichiometry,
        constants=constants,
        dynamic_indices=dynamic_indices,
        repeated_assignment_order=repeated_order,
        initial_assignment_rules=initial_rules,
    )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _snapshot_digest(directory: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(directory.rglob("*")):
        if path.is_dir():
            continue
        digest.update(path.relative_to(directory).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _enumerate_dose_times(dose: DoseEntry, days: float) -> List[float]:
    interval = float(dose.interval)
    repeat = max(int(dose.repeat_count), 0)
    times: List[float] = []
    total = repeat + 1
    for occurrence in range(total):
        time_point = float(dose.start_time) + occurrence * interval
        if time_point > days + _TIME_TOL:
            break
        times.append(time_point)
        if interval <= _TIME_TOL:
            break
    return times


def _fallback_anti_pd1_doses(model: FrozenModel, days: float) -> List[ScheduledDose]:
    mg_per_kg = 3.0
    patient_weight = 70.0
    molecular_weight = 1.436e8
    amount = patient_weight * mg_per_kg / molecular_weight
    fallback = DoseEntry(
        index=10_000,
        name="nivolumab_fallback",
        dose_type="RepeatDose",
        target="V_C.nivolumab",
        amount=amount,
        amount_units="mole",
        start_time=0.0,
        interval=14.0,
        repeat_count=30,
    )
    schedule: List[ScheduledDose] = []
    for time_point in _enumerate_dose_times(fallback, days):
        schedule.append(ScheduledDose(time=time_point, priority=fallback.index, dose=fallback, amount=amount))
    return schedule


def _build_scheduled_doses(
    model: FrozenModel,
    therapy: str,
    days: float,
    *,
    custom_doses: Optional[Sequence[DoseEntry]] = None,
) -> List[ScheduledDose]:
    therapy_flag = therapy.lower()
    scheduled: List[ScheduledDose] = []
    if custom_doses is not None:
        base_doses = list(custom_doses)
    else:
        base_doses = [] if therapy_flag == "none" else list(model.doses)
    for dose in base_doses:
        for time_point in _enumerate_dose_times(dose, days):
            scheduled.append(ScheduledDose(time=time_point, priority=dose.index, dose=dose, amount=dose.amount))
    if not scheduled and therapy_flag == "anti_pd1" and custom_doses is None:
        scheduled.extend(_fallback_anti_pd1_doses(model, days))
    scheduled.sort(key=lambda item: (item.time, item.priority))
    return scheduled


def _record_solution_samples(
    samples: Dict[float, np.ndarray],
    sample_times: np.ndarray,
    start: float,
    end: float,
    sol,
    model: FrozenModel,
    *,
    inclusive: bool,
) -> None:
    if sol.sol is None:
        return
    for t_sample in sample_times:
        if t_sample in samples:
            continue
        if t_sample <= start + _TIME_TOL:
            continue
        if inclusive:
            if t_sample > end + _TIME_TOL:
                continue
        else:
            if t_sample >= end - _TIME_TOL:
                continue
        vec = np.asarray(sol.sol(t_sample), dtype=float)
        ctx = model.build_context_from_state(vec.copy())
        model.evaluate_repeated_assignments(ctx)
        model.apply_algebraic_rules(ctx, vec, mutate=False)
        model.sync_state_from_context(ctx, vec)
        samples[t_sample] = vec


def _perform_t0_quick_check(model: FrozenModel, state: np.ndarray, context: Dict[str, float]) -> None:
    working_state = np.array(state, dtype=float, copy=True)
    working_context = dict(context)
    model.apply_algebraic_rules(working_context, working_state, mutate=False)
    flux_values = model.evaluate_reactions(working_context)
    derivative = model.evaluate_ode_rhs(working_context)
    ode_map = {
        identifier: derivative[idx]
        for identifier, idx in model.dynamic_indices.items()
    }

    evaluation: Dict[str, float] = {}
    for name, value in flux_values.items():
        evaluation[f"flux:{name}"] = float(value)
    for identifier, value in ode_map.items():
        evaluation[f"ode:{identifier}"] = float(value)
    for entry in model.events:
        if entry.trigger_boolean_compiled is None:
            continue
        raw = entry.trigger_boolean_compiled.evaluate_raw(working_context)
        evaluation[f"event:{entry.index}"] = 1.0 if bool(raw) else 0.0

    reference_path = model.source_dir / "equations_eval_t0_reference.csv"
    reference_values: Dict[str, float] = {}
    if reference_path.is_file():
        frame = pd.read_csv(reference_path)
        if "name" in frame.columns:
            value_column = None
            for candidate in ("reference_value", "value", "matlab_value"):
                if candidate in frame.columns:
                    value_column = candidate
                    break
            if value_column is None:
                value_column = "reference_value"
            for row in frame.itertuples(index=False):
                name = str(row.name)
                if not hasattr(row, value_column):
                    continue
                raw_ref = getattr(row, value_column)
                if pd.isna(raw_ref):
                    continue
                reference_values[name] = float(raw_ref)

    rows: List[Dict[str, float]] = []
    failures: List[Tuple[str, float, float, float]] = []
    for name, python_value in sorted(evaluation.items()):
        ref_value = reference_values.get(name)
        python_float = float(python_value)
        if ref_value is None:
            rel_err = float("nan")
        else:
            abs_err = abs(python_float - ref_value)
            rel_err = abs_err / max(abs(ref_value), _T0_ABS_TOL)
            if abs_err > (_T0_ABS_TOL + _T0_REL_TOL * abs(ref_value)):
                failures.append((name, python_float, ref_value, abs_err))
        rows.append(
            {
                "name": name,
                "python_value": python_float,
                "reference_value": ref_value if ref_value is not None else np.nan,
                "rel_err": rel_err,
            }
        )

    output_path = model.source_dir / "equations_eval_t0.csv"
    try:
        pd.DataFrame(rows).to_csv(output_path, index=False)
    except OSError:
        pass

    if reference_values and failures:
        details = ", ".join(f"{name} Δ={abs_err:.2e}" for name, _, _, abs_err in failures[:5])
        raise RuntimeError(f"t=0 quick check failed for snapshot '{model.name}': {details}")


def simulate_frozen_model(
    snapshot: str,
    *,
    days: float,
    therapy: str,
    seed: Optional[int] = None,
    emit_diagnostics: bool = False,
    run_label: Optional[str] = None,
    event_log: Optional[List[Dict[str, object]]] = None,
    rtol_override: Optional[float] = None,
    atol_override: Optional[float] = None,
    sample_interval_hours: Optional[float] = 24.0,
    extra_outputs: Optional[Mapping[str, str]] = None,
    custom_doses: Optional[Sequence[DoseEntry]] = None,
) -> ScenarioResult:
    model = load_frozen_model(snapshot)
    state = model.initial_state().astype(float)
    model.apply_initial_assignments_to_state(state)

    if event_log is not None:
        event_log.clear()

    solver_options = model.config.get("SolverOptions", {})
    rel_tol = solver_options.get("RelativeTolerance", 1e-7)
    abs_tol = solver_options.get("AbsoluteTolerance", 1e-10)
    max_step_value = solver_options.get("MaxStep", None)
    if isinstance(max_step_value, (list, tuple)) and len(max_step_value) == 0:
        max_step_value = None
    max_step = float(max_step_value) if max_step_value not in (None, []) else np.inf
    if max_step == 0:
        max_step = np.inf

    if emit_diagnostics:
        logger = logging.getLogger(__name__)
        time_units = model.config.get("TimeUnits", "day")
        method = model.config.get("SolverType", "BDF")
        max_step_display = "inf" if not math.isfinite(max_step) else f"{max_step:g}"
        equations_path = model.source_dir / "equations.txt"
        config_path = model.source_dir / "configset.json"
        snapshot_sha = _snapshot_digest(model.source_dir)
        equations_sha = _sha256_file(equations_path) if equations_path.exists() else "NA"
        config_sha = _sha256_file(config_path) if config_path.exists() else "NA"
        banner_label = run_label or snapshot
        logger.info(
            "solver=%s rtol=%g atol=%g max_step=%s time_units=%s seed=%s run=%s therapy=%s stop_time=%s",
            method,
            rel_tol,
            abs_tol,
            max_step_display,
            time_units,
            seed,
            banner_label,
            therapy,
            days,
        )
        logger.info(
            "semantics=%s snapshot_sha=%s equations_sha=%s configset_sha=%s",
            SEMANTICS_VERSION,
            snapshot_sha,
            equations_sha,
            config_sha,
        )

    if rtol_override is not None:
        rel_tol = rtol_override
    if atol_override is not None:
        abs_tol = atol_override

    if sample_interval_hours is not None and sample_interval_hours <= 0:
        raise ValueError("sample_interval_hours must be positive")

    def _build_sample_times(stop_time: float, interval_hours: Optional[float]) -> np.ndarray:
        tol = 1e-12
        if interval_hours is None:
            last_complete_day = int(math.floor(stop_time + tol))
            times = np.arange(0, last_complete_day + 1, dtype=float)
            if abs(stop_time - float(last_complete_day)) > tol:
                times = np.append(times, stop_time)
            elif times[-1] != stop_time:
                times = np.append(times, stop_time)
            return np.unique(np.round(times, 12))
        step_days = float(interval_hours) / 24.0
        if step_days <= tol:
            step_days = tol
        max_steps = int(math.floor((stop_time + tol) / step_days))
        times = np.arange(0.0, (max_steps + 1) * step_days + tol, step_days)
        times = np.round(times, 12)
        if not np.isclose(times[0], 0.0, atol=tol):
            times = np.concatenate(([0.0], times))
        if times[-1] < stop_time - tol:
            times = np.append(times, round(stop_time, 12))
        else:
            times[-1] = round(stop_time, 12)
        return np.unique(times)

    sample_times = _build_sample_times(days, sample_interval_hours)
    tol_time = _TIME_TOL

    def reconcile(vec: np.ndarray) -> Dict[str, float]:
        ctx = model.build_context_from_state(vec.copy())
        model.evaluate_repeated_assignments(ctx)
        model.apply_algebraic_rules(ctx, vec, mutate=False)
        model.sync_state_from_context(ctx, vec)
        return ctx

    def record_sample(samples: Dict[float, np.ndarray], time: float, vec: np.ndarray) -> None:
        for target in sample_times:
            if abs(time - target) <= tol_time:
                samples[target] = vec.copy()
                break

    context = reconcile(state)
    _perform_t0_quick_check(model, state, context)

    samples: Dict[float, np.ndarray] = {}
    record_sample(samples, 0.0, state)

    scheduled_doses = _build_scheduled_doses(model, therapy, days, custom_doses=custom_doses)
    dose_index = 0
    pending_events: List[ScheduledEvent] = []

    def make_event_function(entry: EventEntry):
        def fn(t, y):
            vec = np.asarray(y, dtype=float)
            ctx = model.build_context_from_state(vec.copy())
            model.evaluate_repeated_assignments(ctx)
            model.apply_algebraic_rules(ctx, vec, mutate=False)
            return entry.trigger_compiled.evaluate(ctx)

        fn.direction = entry.direction
        fn.terminal = False
        return fn

    event_functions = [make_event_function(entry) for entry in model.events]

    current_state = state
    context = reconcile(current_state)

    while dose_index < len(scheduled_doses) and abs(scheduled_doses[dose_index].time) <= tol_time:
        scheduled_dose = scheduled_doses[dose_index]
        model.apply_dose(scheduled_dose.dose, scheduled_dose.amount, context, current_state)
        context = reconcile(current_state)
        record_sample(samples, 0.0, current_state)
        dose_index += 1

    current_time = 0.0

    while current_time < days - tol_time:
        next_dose_time = scheduled_doses[dose_index].time if dose_index < len(scheduled_doses) else float("inf")
        next_event_time = pending_events[0].time if pending_events else float("inf")
        target_time = min(days, next_dose_time, next_event_time)

        if target_time <= current_time + tol_time:
            # Process any queued actions at this instant without integrating.
            current_time = target_time
            while pending_events and abs(pending_events[0].time - current_time) <= tol_time:
                scheduled_event = heapq.heappop(pending_events)
                for assignment in scheduled_event.entry.assignments:
                    value = assignment.compiled.evaluate(context) if assignment.compiled is not None else 0.0
                    model._apply_target_value(assignment.target, value, context, current_state)
                context = reconcile(current_state)
                record_sample(samples, current_time, current_state)
            while dose_index < len(scheduled_doses) and scheduled_doses[dose_index].time <= current_time + tol_time:
                scheduled_dose = scheduled_doses[dose_index]
                model.apply_dose(scheduled_dose.dose, scheduled_dose.amount, context, current_state)
                context = reconcile(current_state)
                record_sample(samples, current_time, current_state)
                dose_index += 1
            continue

        sol = solve_ivp(
            model.rhs,
            (current_time, target_time),
            current_state,
            method="BDF",
            rtol=rel_tol,
            atol=abs_tol,
            max_step=max_step,
            dense_output=True,
            events=event_functions if event_functions else None,
            vectorized=False,
        )

        if not sol.success:  # pragma: no cover - propagates solver failure
            raise RuntimeError(f"Integration failed at t={current_time}: {sol.message}")

        triggered_pairs: List[Tuple[float, EventEntry]] = []
        if sol.t_events and event_functions:
            for idx, times in enumerate(sol.t_events):
                for time_val in np.atleast_1d(times):
                    triggered_pairs.append((float(time_val), model.events[idx]))

        if triggered_pairs:
            triggered_pairs.sort(key=lambda pair: (pair[0], pair[1].index))
            event_time = triggered_pairs[0][0]
            _record_solution_samples(samples, sample_times, current_time, event_time, sol, model, inclusive=False)
            if sol.sol is not None:
                current_state = np.asarray(sol.sol(event_time), dtype=float)
            else:
                idx = int(np.argmin(np.abs(sol.t - event_time)))
                current_state = np.asarray(sol.y[:, idx], dtype=float)
            context = reconcile(current_state)
            same_time_entries = [entry for time_val, entry in triggered_pairs if abs(time_val - event_time) <= tol_time]
            for entry in same_time_entries:
                assignments_text = entry.assignments_text
                delay_value = 0.0
                if entry.delay_compiled is not None:
                    delay_value = float(entry.delay_compiled.evaluate(context))
                if entry.delay_type == "time" and delay_value > tol_time:
                    scheduled_time = event_time + delay_value
                    heapq.heappush(
                        pending_events,
                        ScheduledEvent(
                            time=scheduled_time,
                            priority=entry.index,
                            entry=entry,
                            trigger_time=event_time,
                            delay=delay_value,
                            assignments=assignments_text,
                        ),
                    )
                    continue
                if event_log is not None:
                    event_log.append(
                        {
                            "event_index": entry.index,
                            "time_fire": event_time,
                            "time_trigger": event_time,
                            "delay": delay_value,
                            "type": "immediate",
                            "assignments": assignments_text,
                        }
                    )
                for assignment in entry.assignments:
                    value = assignment.compiled.evaluate(context) if assignment.compiled is not None else 0.0
                    model._apply_target_value(assignment.target, value, context, current_state)
                context = reconcile(current_state)
            record_sample(samples, event_time, current_state)
            current_time = event_time
        else:
            t_end = float(sol.t[-1])
            _record_solution_samples(samples, sample_times, current_time, t_end, sol, model, inclusive=True)
            current_state = np.asarray(sol.y[:, -1], dtype=float)
            context = reconcile(current_state)
            current_time = t_end

        while pending_events and abs(pending_events[0].time - current_time) <= tol_time:
            scheduled_event = heapq.heappop(pending_events)
            if event_log is not None:
                event_log.append(
                    {
                        "event_index": scheduled_event.entry.index,
                        "time_fire": current_time,
                        "time_trigger": scheduled_event.trigger_time,
                        "delay": scheduled_event.delay,
                        "type": "delayed",
                        "assignments": scheduled_event.assignments,
                    }
                )
            for assignment in scheduled_event.entry.assignments:
                value = assignment.compiled.evaluate(context) if assignment.compiled is not None else 0.0
                model._apply_target_value(assignment.target, value, context, current_state)
            context = reconcile(current_state)
            record_sample(samples, current_time, current_state)

        while dose_index < len(scheduled_doses) and scheduled_doses[dose_index].time <= current_time + tol_time:
            scheduled_dose = scheduled_doses[dose_index]
            model.apply_dose(scheduled_dose.dose, scheduled_dose.amount, context, current_state)
            context = reconcile(current_state)
            record_sample(samples, current_time, current_state)
            dose_index += 1

    record_sample(samples, current_time, current_state)

    output_states = []
    last_state = current_state.copy()
    for t_val in sample_times:
        if t_val in samples:
            last_state = samples[t_val]
        output_states.append(last_state.copy())
    states = np.vstack(output_states)

    cancer_cells = []
    dead_cells = []
    t_cells = []
    tumour_volume_l = []
    tumour_diameter_cm = []
    pd1_occupancy = []
    tcell_density = []
    extra_columns = extra_outputs or {}
    extras_accumulators: Dict[str, List[float]] = {name: [] for name in extra_columns}

    for time_point, vector in zip(sample_times, states):
        vec = vector.copy()
        context = model.build_context_from_state(vec)
        model.evaluate_repeated_assignments(context)
        model.apply_algebraic_rules(context, vec, mutate=False)
        model.sync_state_from_context(context, vec)

        c_cells = context.get("C1", 0.0)
        d_cells = context.get("C_x", 0.0)
        t_tumour = context.get("V_T.T1", 0.0)
        tumour_t0 = context.get("V_T.T0", 0.0)
        volume_l = context.get("V_T", 0.0)
        diameter_cm = ((3.0 * volume_l * 1e3) / (4.0 * math.pi)) ** (1.0 / 3.0) * 2.0 if volume_l > 0 else 0.0
        occupancy = context.get("H_PD1_C1", 0.0)
        density = t_tumour / max(volume_l * 1e6, 1e-12)
        t_total = context.get("T_total", t_tumour + tumour_t0)

        cancer_cells.append(c_cells)
        dead_cells.append(d_cells)
        t_cells.append(t_total)
        tumour_volume_l.append(volume_l)
        tumour_diameter_cm.append(diameter_cm)
        pd1_occupancy.append(occupancy)
        tcell_density.append(density)
        for column, context_key in extra_columns.items():
            extras_accumulators[column].append(float(context.get(context_key, 0.0)))

    extras_arrays = {
        column: np.array(values, dtype=float) for column, values in extras_accumulators.items()
    }

    return ScenarioResult(
        time_days=sample_times,
        cancer_cells=np.array(cancer_cells, dtype=float),
        dead_cells=np.array(dead_cells, dtype=float),
        t_cells=np.array(t_cells, dtype=float),
        tumour_volume_l=np.array(tumour_volume_l, dtype=float),
        tumour_diameter_cm=np.array(tumour_diameter_cm, dtype=float),
        pd1_occupancy=np.array(pd1_occupancy, dtype=float),
        tcell_density_per_ul=np.array(tcell_density, dtype=float),
        extras=extras_arrays,
    )
