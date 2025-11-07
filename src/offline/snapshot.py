"""Snapshot loading and model construction for frozen SimBiology projects."""

from __future__ import annotations

import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import sympy as sp

from .entities import (
    CompiledExpression,
    DoseEntry,
    EventAssignment,
    EventEntry,
    OdeEntry,
    ReactionEntry,
    RuleEntry,
    SpeciesEntry,
)
from .errors import SnapshotError
from .units import convert_amount, convert_parameter_value, convert_volume


class FrozenModel:
    """Immutable representation of a frozen SimBiology snapshot."""

    def __init__(
        self,
        *,
        name: str,
        source_dir: Path,
        species: List[SpeciesEntry],
        species_lookup: Dict[str, SpeciesEntry],
        species_name_lookup: Dict[str, SpeciesEntry],
        parameters: Dict[str, float],
        compartments: Dict[str, float],
        config: Dict[str, object],
        rules: List[RuleEntry],
        rate_rules: List[RuleEntry],
        algebraic_rules: List[RuleEntry],
        reactions: List[ReactionEntry],
        odes: List[OdeEntry],
        events: List[EventEntry],
        doses: List[DoseEntry],
        variants: List[Dict[str, object]],
        stoichiometry: List[List[Tuple[int, float]]],
        constants: Dict[str, float],
        dynamic_indices: Dict[str, int],
        repeated_assignment_order: List[RuleEntry],
        initial_assignment_rules: List[RuleEntry],
        time_unit: str,
        solver_type: str,
        provenance: Dict[str, str],
    ) -> None:
        self.name = name
        self.source_dir = source_dir
        self.species = species
        self.species_lookup = species_lookup
        self.species_name_lookup = species_name_lookup
        self.parameters = parameters
        self.compartments = compartments
        self.config = config
        self.rules = rules
        self.rate_rules = rate_rules
        self.algebraic_rules = algebraic_rules
        self.reactions = reactions
        self.odes = odes
        self.events = events
        self.doses = doses
        self.variants = variants
        self.stoichiometry = stoichiometry
        self.constants = constants
        self.dynamic_indices = dynamic_indices
        self.repeated_assignment_order = repeated_assignment_order
        self.initial_assignment_rules = initial_assignment_rules
        self.time_unit = time_unit
        self.solver_type = solver_type
        self.provenance = provenance
        self._jacobian_sparsity: Optional[np.ndarray] = None

    # --- state/context manipulation -------------------------------------------------

    def initial_state(self) -> np.ndarray:
        state = np.zeros(len(self.dynamic_indices), dtype=float)
        for species in self.species:
            if species.identifier in self.dynamic_indices:
                state[self.dynamic_indices[species.identifier]] = species.initial_amount
        return state

    def _base_context(self) -> Dict[str, float]:
        ctx: Dict[str, float] = {}
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

    def _resolve_dynamic_index(self, symbol: str) -> Optional[int]:
        idx = self.dynamic_indices.get(symbol)
        if idx is not None:
            return idx
        entry = self.species_lookup.get(symbol)
        if entry and entry.identifier in self.dynamic_indices:
            return self.dynamic_indices[entry.identifier]
        named = self.species_name_lookup.get(symbol)
        if named and named.identifier in self.dynamic_indices:
            return self.dynamic_indices[named.identifier]
        return None

    def _build_jacobian_sparsity(self) -> np.ndarray:
        size = len(self.dynamic_indices)
        pattern = np.zeros((size, size), dtype=bool)

        def mark(row_identifier: str, tokens: Sequence[str]) -> None:
            row_idx = self._resolve_dynamic_index(row_identifier)
            if row_idx is None:
                return
            for token in tokens:
                col_idx = self._resolve_dynamic_index(token)
                if col_idx is not None:
                    pattern[row_idx, col_idx] = True

        for ode in self.odes:
            if ode.compiled is None:
                continue
            mark(ode.identifier, ode.compiled.tokens)

        for rule in self.rate_rules:
            if rule.compiled is None:
                continue
            mark(rule.target, rule.compiled.tokens)

        np.fill_diagonal(pattern, True)
        return pattern

    def jacobian_sparsity(self) -> np.ndarray:
        if self._jacobian_sparsity is None:
            self._jacobian_sparsity = self._build_jacobian_sparsity()
        return self._jacobian_sparsity

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

    # --- evaluation helpers --------------------------------------------------------

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
            entry = self.species_lookup.get(ode.identifier)
            if target_idx is None:
                name_entry = self.species_name_lookup.get(ode.identifier)
                if name_entry is not None:
                    entry = name_entry
                    target_idx = self.dynamic_indices.get(name_entry.identifier)
            if target_idx is None:
                continue
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

    def apply_dose(self, dose: DoseEntry, amount: float, context: Dict[str, float], state: np.ndarray) -> Dict[str, object]:
        target = dose.target
        amount_units = dose.amount_units or "mole"
        mw = getattr(dose, "molecular_weight_g_per_mol", None)
        try:
            amount_mol = convert_amount(amount, amount_units, mw)
        except ValueError:
            amount_mol = amount
        entry = self.species_lookup.get(target)
        if entry is None and target in self.species_name_lookup:
            entry = self.species_name_lookup[target]
        if entry is None:
            new_value = context.get(target, 0.0) + amount_mol
            self._apply_target_value(target, new_value, context, state)
            return {
                "target": target,
                "interpreted_dimension": "amount",
                "units": "mole",
                "compartment": None,
                "compartment_volume_l": None,
                "delta_state_value": amount_mol,
                "delta_amount_mol": amount_mol,
            }
        delta = amount_mol
        units_lower = (entry.units or "").lower()
        dimension_lower = (entry.interpreted_dimension or "").lower()

        def _looks_like_concentration() -> bool:
            if dimension_lower:
                if any(token in dimension_lower for token in ("concentration", "amount/vol", "mass/vol", "mol/vol")):
                    return True
                if dimension_lower in {"amount", "substance", "mass"}:
                    return False
            if not units_lower:
                return False
            if "molar" in units_lower or "mol/" in units_lower or "mole/" in units_lower:
                return True
            if "/" in units_lower:
                return True
            return False

        is_concentration = _looks_like_concentration()
        if is_concentration:
            compartment_name = entry.compartment
            compartment_volume = context.get(compartment_name, self.compartments.get(compartment_name))
            if compartment_volume in (None, 0.0):
                raise RuntimeError(f"Missing compartment volume for dose target {target}")
            delta = amount / compartment_volume
        def _label(value: Optional[str], fallback: str) -> str:
            if value is None:
                return fallback
            if isinstance(value, float) and math.isnan(value):
                return fallback
            text = str(value).strip()
            return text or fallback

        dimension_label = _label(entry.interpreted_dimension, "molarity" if is_concentration else "amount")
        units_label = _label(entry.units, "molarity" if is_concentration else "mole")
        current_value = context.get(entry.identifier, 0.0)
        new_value = current_value + delta
        self._apply_target_value(entry.identifier, new_value, context, state)
        return {
            "target": entry.identifier,
            "interpreted_dimension": dimension_label,
            "units": units_label,
            "compartment": entry.compartment,
            "compartment_volume_l": context.get(entry.compartment, self.compartments.get(entry.compartment)),
            "delta_state_value": delta,
            "delta_amount_mol": amount_mol,
        }

    def rhs(self, t: float, y: np.ndarray) -> np.ndarray:
        state_view = np.array(y, dtype=float, copy=True)
        context = self.build_context_from_state(state_view)
        self.evaluate_repeated_assignments(context)
        self.apply_algebraic_rules(context, state_view, mutate=False)
        self.evaluate_reactions(context)
        derivative = self.evaluate_ode_rhs(context)
        return derivative


# --------------------------------------------------------------------------------------
# Snapshot loading helpers
# --------------------------------------------------------------------------------------


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


def _load_rows(path: Path) -> List[Dict[str, object]]:
    suffix = path.suffix.lower()
    if suffix == ".json":
        with path.open("r", encoding="utf8") as handle:
            data = json.load(handle)
        if isinstance(data, list):
            return data
        raise SnapshotError(f"Unsupported JSON structure in {path}")
    frame = _read_csv(path)
    return [row._asdict() for row in frame.itertuples(index=False)]


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
    rows = _load_rows(path)
    base_params: Dict[str, float] = {}
    derived_specs: List[Dict[str, object]] = []
    order: List[str] = []
    for row in rows:
        name = str(row.get("name"))
        units = str(row.get("units", "")) if row.get("units") is not None else ""
        value = row.get("value")
        if value is None:
            derived_specs.append(row)
            continue
        converted = convert_parameter_value(float(value), units)
        base_params[name] = converted
        order.append(name)
    resolved = _evaluate_derived_parameters(base_params, derived_specs)
    for spec in derived_specs:
        name = str(spec.get("name"))
        if name in resolved and name not in order:
            order.append(name)
    return resolved, order


def _load_compartments(path: Path) -> Tuple[Dict[str, float], List[str]]:
    rows = _load_rows(path)
    comps: Dict[str, float] = {}
    order: List[str] = []
    for row in rows:
        name = str(row.get("name"))
        value = float(row.get("capacity", 0.0))
        units = str(row.get("units", "")) if row.get("units") is not None else ""
        comps[name] = _convert_compartment_value(value, units)
        order.append(name)
    return comps, order


def _convert_compartment_value(value: float, units: str) -> float:
    return convert_volume(value, units)


def _convert_parameter_value(value: float, units: str) -> float:
    return convert_parameter_value(value, units)


def _evaluate_derived_parameters(
    base_params: Dict[str, float],
    derived_specs: Sequence[Dict[str, object]],
) -> Dict[str, float]:
    resolved = dict(base_params)
    remaining = list(derived_specs)

    def _evaluate_expression(expr: str, values: Sequence[float]) -> float:
        substituted = expr
        for idx, val in enumerate(values, start=1):
            substituted = substituted.replace(f"p({idx})", f"({val})")
        safe_locals = {"__builtins__": None, "pi": math.pi, "e": math.e}
        return float(eval(substituted, safe_locals, {}))

    while remaining:
        progress = False
        for spec in list(remaining):
            name = str(spec.get("name"))
            deps = [str(dep) for dep in spec.get("derived_from", []) or []]
            if any(dep not in resolved for dep in deps):
                continue
            expr = spec.get("expression")
            units = str(spec.get("units", "")) if spec.get("units") is not None else ""
            if expr:
                values = [resolved[dep] for dep in deps]
                raw_value = _evaluate_expression(expr, values)
            elif spec.get("value") is not None:
                raw_value = float(spec.get("value"))
            else:
                continue
            resolved[name] = convert_parameter_value(raw_value, units)
            remaining.remove(spec)
            progress = True
        if not progress:
            unresolved = [str(spec.get("name")) for spec in remaining]
            raise SnapshotError(f"Unable to resolve derived parameters: {unresolved}")
    return resolved


def _load_rules(
    path: Path,
    tokens_sorted: Sequence[str],
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
                compiled = _compile_expression(raw_expr, tokens_sorted, safe_names, reverse_safe)
        dep_tokens: Tuple[str, ...] = tuple()
        if compiled is not None:
            dep_tokens = tuple(token for token in compiled.tokens if token != target)
        safe_target = safe_names.get(target) if target else None
        algebraic_solver: Optional[CompiledExpression] = None
        if rule_type_lower in {"algebraic", "algebraicrule"} and compiled is not None and safe_target:
            try:
                target_symbol = sp.Symbol(safe_target)
                solutions = sp.solve(sp.Eq(compiled.sympy_expr, 0), target_symbol, dict=True)
            except Exception:  # pragma: no cover
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


def _load_reactions(
    path: Path,
    tokens_sorted: Sequence[str],
    safe_names: Dict[str, str],
    reverse_safe: Dict[str, str],
    fluxes: Dict[str, str],
) -> List[ReactionEntry]:
    frame = _read_csv(path)
    reactions: List[ReactionEntry] = []
    for row in frame.itertuples(index=False):
        expression = str(getattr(row, "kinetic_expression", ""))
        if not expression or expression.lower() == "nan":
            expression = str(getattr(row, "reaction_rate", ""))
        reaction_name = str(row.name)
        if reaction_name in fluxes:
            expression = fluxes[reaction_name]
        elif expression in fluxes:
            expression = fluxes[expression]
        compiled = _compile_expression(expression, tokens_sorted, safe_names, reverse_safe)
        reactions.append(
            ReactionEntry(
                name=str(row.name),
                expression=expression,
                compiled=compiled,
            )
        )
    return reactions


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


def _load_events(
    path: Path,
    tokens_sorted: Sequence[str],
    safe_names: Dict[str, str],
    reverse_safe: Dict[str, str],
) -> List[EventEntry]:
    frame = _read_csv(path)
    events: List[EventEntry] = []
    for row in frame.itertuples(index=False):
        trigger_expr_raw = str(row.trigger)
        expr_str, direction = _parse_trigger(trigger_expr_raw)
        trigger_compiled = _compile_expression(expr_str, tokens_sorted, safe_names, reverse_safe)
        trigger_bool_compiled: Optional[CompiledExpression] = None
        if trigger_expr_raw:
            trigger_bool_compiled = _compile_expression(trigger_expr_raw, tokens_sorted, safe_names, reverse_safe)
        delay_expr = "" if pd.isna(row.delay) else str(row.delay)
        delay_compiled = _compile_expression(delay_expr, tokens_sorted, safe_names, reverse_safe) if delay_expr else None
        delay_type_raw = "" if not hasattr(row, "delay_type") or pd.isna(row.delay_type) else str(row.delay_type)
        delay_type = delay_type_raw.strip().lower() or "time"
        assignments_raw = _parse_assignments(str(row.assignments))
        assignments: List[EventAssignment] = []
        for target, expr in assignments_raw:
            if expr:
                compiled = _compile_expression(expr, tokens_sorted, safe_names, reverse_safe)
            else:
                compiled = None
            assignments.append(EventAssignment(target=target, expression=expr, compiled=compiled))
        events.append(
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
    return events


def _load_doses(path: Path) -> List[DoseEntry]:
    rows = _load_rows(path)
    doses: List[DoseEntry] = []
    for row in rows:
        amount = float(row.get("amount", 0.0) or 0.0)
        amount_units = str(row.get("amount_units", "")) if row.get("amount_units") is not None else ""
        rate = row.get("rate")
        rate_units = str(row.get("rate_units", "")) if row.get("rate_units") is not None else ""
        duration = row.get("duration")
        amount_mg = row.get("amount_mg")
        mw = row.get("molecular_weight_g_per_mol")
        doses.append(
            DoseEntry(
                index=int(row.get("dose_index", len(doses) + 1)),
                name=str(row.get("name")),
                dose_type=str(row.get("type", row.get("dose_type", ""))),
                target=str(row.get("target")),
                amount=amount,
                amount_units=amount_units,
                start_time=float(row.get("start_time", 0.0)),
                interval=float(row.get("interval", 0.0)),
                repeat_count=int(row.get("repeat_count", 0)),
                rate=float(rate) if rate not in (None, "") else None,
                rate_units=rate_units,
                duration=float(duration) if duration not in (None, "") else None,
                amount_mg=float(amount_mg) if amount_mg not in (None, "") else None,
                molecular_weight_g_per_mol=float(mw) if mw not in (None, "") else None,
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
                "content": str(getattr(row, "content", "")),
            }
        )
    return variants


def _load_stoichiometry(
    path: Path,
    reaction_names: Sequence[str],
    dynamic_indices: Dict[str, int],
) -> List[List[Tuple[int, float]]]:
    frame = _read_csv(path)
    name_to_idx = {name: idx for idx, name in enumerate(reaction_names)}
    stoich: List[List[Tuple[int, float]]] = [[] for _ in reaction_names]
    for row in frame.itertuples(index=False):
        reaction_idx = name_to_idx.get(str(row.reaction))
        species_idx = dynamic_indices.get(str(row.species))
        if reaction_idx is None or species_idx is None:
            continue
        coeff = getattr(row, "coefficient", None)
        if coeff is None or (isinstance(coeff, float) and math.isnan(coeff)):
            coeff = getattr(row, "stoichiometry", 0.0)
        stoich[reaction_idx].append((species_idx, float(coeff)))
    return stoich


def _topological_order_repeated(rules: List[RuleEntry]) -> List[RuleEntry]:
    dependencies: Dict[str, List[str]] = {}
    for rule in rules:
        dependencies[rule.target] = list(rule.dependencies)
    ordered: List[RuleEntry] = []
    temporary: Dict[str, bool] = {}
    permanent: Dict[str, bool] = {}
    rule_by_target = {rule.target: rule for rule in rules}

    def visit(node: str) -> None:
        if permanent.get(node):
            return
        if temporary.get(node):
            raise RuntimeError(f"Cycle detected in repeated assignments at '{node}'")
        temporary[node] = True
        for dependency in dependencies.get(node, []):
            visit(dependency)
        temporary.pop(node, None)
        permanent[node] = True
        if node in rule_by_target:
            ordered.append(rule_by_target[node])

    for rule in rules:
        visit(rule.target)
    return ordered


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot_digest(directory: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(directory.rglob("*")):
        if path.is_dir():
            continue
        digest.update(path.relative_to(directory).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


SnapshotSpecifier = Union[str, Path]

REQUIRED_SNAPSHOT_FILES = [
    "configset.json",
    "equations.txt",
    "species.csv",
    "parameters.csv",
    "compartments.csv",
    "rules.csv",
    "reactions.csv",
    "events.csv",
    "doses.csv",
    "stoichiometry.csv",
]


def _resolve_snapshot_path(specifier: SnapshotSpecifier, root: Optional[Path]) -> Path:
    if isinstance(specifier, Path):
        candidate = specifier
    else:
        candidate = Path(specifier)
    if candidate.is_dir():
        return candidate
    if candidate.exists() and not candidate.is_dir():
        raise SnapshotError(f"Snapshot path '{candidate}' exists but is not a directory")
    base_root = root if root is not None else Path("artifacts") / "matlab_frozen_model"
    resolved = base_root / str(specifier)
    if resolved.is_dir():
        return resolved
    raise SnapshotError(f"Snapshot directory {resolved} is missing. Run export_matlab_snapshot first.")


def load_frozen_model(specifier: SnapshotSpecifier, *, root: Optional[Path] = None) -> FrozenModel:
    model_dir = _resolve_snapshot_path(specifier, root).resolve()
    return _load_frozen_model_cached(str(model_dir))


@lru_cache(maxsize=None)
def _load_frozen_model_cached(model_dir_str: str) -> FrozenModel:
    model_dir = Path(model_dir_str)
    name = model_dir.name
    missing = [entry for entry in REQUIRED_SNAPSHOT_FILES if not (model_dir / entry).is_file()]
    if missing:
        raise SnapshotError(f"Snapshot '{model_dir}' is missing required files: {', '.join(missing)}")
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

    equations_path = model_dir / "equations.txt"
    config_path = model_dir / "configset.json"
    parameters_path = model_dir / "parameters.csv"

    fluxes = _load_flux_expressions(equations_path)

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

    time_unit = str(config.get("TimeUnits", "") or "")
    solver_type = str(config.get("SolverType", "BDF") or "BDF")
    species_units = {entry.identifier: entry.units for entry in species if entry.units}
    provenance = {
        "snapshot_name": name,
        "snapshot_sha": snapshot_digest(model_dir),
        "equations_sha": sha256_file(equations_path) if equations_path.exists() else "",
        "config_sha": sha256_file(config_path) if config_path.exists() else "",
        "parameters_sha": sha256_file(parameters_path) if parameters_path.exists() else "",
        "time_unit": time_unit,
        "solver_type": solver_type,
        "sympy_version": sp.__version__,
        "parser_flags": "default_v1",
        "units_map": json.dumps(species_units, sort_keys=True),
    }

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
        time_unit=time_unit,
        solver_type=solver_type,
        provenance=provenance,
    )
