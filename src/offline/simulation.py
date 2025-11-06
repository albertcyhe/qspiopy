"""Numerical runtime for frozen SimBiology snapshots."""

from __future__ import annotations

import hashlib
import heapq
import json
import logging
import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Protocol

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .entities import BASE_HEADER, DoseEntry, EventEntry, ScenarioResult, ScheduledDose, SEMANTICS_VERSION
from .errors import AlignmentFail, ConfigError, NumericsError, SnapshotError
from .snapshot import FrozenModel, load_frozen_model

logger = logging.getLogger(__name__)

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


class ExtraOutputs(Protocol):
    """Callable extension point for computing derived outputs."""

    name: str

    def compute(
        self,
        time_days: np.ndarray,
        states: np.ndarray,
        contexts: Sequence[Mapping[str, float]],
        meta: Mapping[str, object],
    ) -> Mapping[str, np.ndarray]:
        ...


class DoseScheduler(Protocol):
    """Strategy for generating scheduled doses."""

    def build(
        self,
        base_doses: List[DoseEntry],
        *,
        model: FrozenModel,
        days: float,
    ) -> Iterable[ScheduledDose]:
        ...


@dataclass(frozen=True)
class SolverConfig:
    """Configuration driving scipy's solve_ivp."""

    method: str
    rtol: float
    atol: float
    max_step: float
    seed: Optional[int]

    def as_dict(self) -> Dict[str, object]:
        return {
            "method": self.method,
            "rtol": self.rtol,
            "atol": self.atol,
            "max_step": self.max_step,
            "seed": self.seed,
        }

    def identity(self) -> str:
        payload = json.dumps(self.as_dict(), sort_keys=True)
        return hashlib.sha256(payload.encode("utf8")).hexdigest()


@dataclass
class ContextKeyOutputs:
    """Adapter for legacy context-key extraction."""

    mapping: Mapping[str, str]

    @property
    def name(self) -> str:
        return "context_keys"

    def compute(
        self,
        time_days: np.ndarray,
        states: np.ndarray,
        contexts: Sequence[Mapping[str, float]],
        meta: Mapping[str, object],
    ) -> Mapping[str, np.ndarray]:
        results: Dict[str, np.ndarray] = {}
        for column, key in self.mapping.items():
            results[column] = np.array([ctx.get(key, 0.0) for ctx in contexts], dtype=float)
        return results


@dataclass(frozen=True, order=True)
class ScheduledEvent:
    time: float
    priority: int
    entry: EventEntry = field(compare=False)
    trigger_time: float = field(compare=False)
    delay: float = field(compare=False)
    assignments: str = field(compare=False)


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
    amount_mg = patient_weight * mg_per_kg
    amount = amount_mg / molecular_weight
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
        amount_mg=amount_mg,
    )
    schedule: List[ScheduledDose] = []
    for time_point in _enumerate_dose_times(fallback, days):
        schedule.append(
            ScheduledDose(
                time=time_point,
                priority=fallback.index,
                dose=fallback,
                amount=amount,
                amount_mg=amount_mg,
            )
        )
    return schedule


def _materialise_schedule(doses: Sequence[DoseEntry], days: float) -> List[ScheduledDose]:
    schedule: List[ScheduledDose] = []
    for dose in doses:
        for time_point in _enumerate_dose_times(dose, days):
            schedule.append(
                ScheduledDose(
                    time=time_point,
                    priority=dose.index,
                    dose=dose,
                    amount=dose.amount,
                    amount_mg=dose.amount_mg,
                )
            )
    schedule.sort(key=lambda item: (item.time, item.priority))
    return schedule


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

    failures: List[Tuple[str, float, float, float]] = []
    for name, python_value in sorted(evaluation.items()):
        ref_value = reference_values.get(name)
        python_float = float(python_value)
        if ref_value is None:
            continue
        abs_err = abs(python_float - ref_value)
        if abs_err > (_T0_ABS_TOL + _T0_REL_TOL * abs(ref_value)):
            failures.append((name, python_float, ref_value, abs_err))

    if failures:
        details = ", ".join(f"{name} Δ={abs_err:.2e}" for name, _, _, abs_err in failures[:5])
        raise AlignmentFail(f"t=0 quick check failed for snapshot '{model.name}': {details}")


def _sample_times(days: float, interval_hours: Optional[float]) -> np.ndarray:
    tol = 1e-12
    if interval_hours is None:
        last_complete_day = int(math.floor(days + tol))
        times = np.arange(0, last_complete_day + 1, dtype=float)
        if abs(days - float(last_complete_day)) > tol:
            times = np.append(times, days)
        elif times[-1] != days:
            times = np.append(times, days)
        return np.unique(np.round(times, 12))
    step_days = float(interval_hours) / 24.0
    if step_days <= tol:
        step_days = tol
    max_steps = int(math.floor((days + tol) / step_days))
    times = np.arange(0.0, (max_steps + 1) * step_days + tol, step_days)
    times = np.round(times, 12)
    if not np.isclose(times[0], 0.0, atol=tol):
        times = np.concatenate(([0.0], times))
    if times[-1] < days - tol:
        times = np.append(times, round(days, 12))
    else:
        times[-1] = round(days, 12)
    return np.unique(times)


def _log_solver_banner(model: FrozenModel, solver: SolverConfig, days: float, emit: bool) -> None:
    if not emit:
        return
    meta = {
        "solver": solver.method,
        "rtol": solver.rtol,
        "atol": solver.atol,
        "max_step": solver.max_step,
        "seed": solver.seed,
        "time_unit": model.time_unit,
        "stop_time": days,
    }
    logger.info("solver_config %s", json.dumps(meta, sort_keys=True))
    logger.info("provenance %s", json.dumps(model.provenance, sort_keys=True))


def _time_unit_sniff(model: FrozenModel, interval_hours: Optional[float]) -> None:
    if interval_hours is None:
        return
    unit = model.time_unit.lower()
    if unit.startswith("day"):
        if abs(interval_hours - 1.0) < 1e-9:
            logger.warning("sampling_interval hours=1 while model time unit is days; verify unit alignment")
        if abs(interval_hours - 24.0) < 1e-9:
            logger.warning("sampling_interval hours=24 while model time unit is days; test for hour/day confusion")
    if unit.startswith("hour"):
        if abs(interval_hours - 24.0) < 1e-9 or abs(interval_hours - (1.0 / 24.0)) < 1e-9:
            logger.warning("hour-based model with interval %.6g may indicate day/hour mismatch", interval_hours)


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
    extra_outputs: Optional[Sequence[ExtraOutputs]] = None,
    scheduler: Optional[DoseScheduler] = None,
    custom_doses: Optional[Sequence[DoseEntry]] = None,
    context_outputs: Optional[Mapping[str, str]] = None,
    dose_audit: Optional[List[Dict[str, object]]] = None,
) -> ScenarioResult:
    try:
        model = load_frozen_model(snapshot)
    except FileNotFoundError as exc:  # pragma: no cover - passthrough
        raise SnapshotError(str(exc)) from exc
    except Exception as exc:  # pragma: no cover
        raise SnapshotError(f"Failed to load snapshot '{snapshot}': {exc}") from exc

    solver_options = model.config.get("SolverOptions", {}) or {}
    rel_tol = rtol_override if rtol_override is not None else solver_options.get("RelativeTolerance", 1e-7)
    abs_tol = atol_override if atol_override is not None else solver_options.get("AbsoluteTolerance", 1e-10)
    max_step_value = solver_options.get("MaxStep", None)
    if isinstance(max_step_value, (list, tuple)) and len(max_step_value) == 0:
        max_step_value = None
    max_step = float(max_step_value) if max_step_value not in (None, []) else np.inf
    if max_step == 0:
        max_step = np.inf

    def _map_solver_type(label: str) -> str:
        lower = label.lower()
        if lower in {"ode15s", "ode23s"}:
            return "BDF"
        if lower == "ode45":
            return "RK45"
        if lower == "ode113":
            return "DOP853"
        return label.upper()

    solver_config = SolverConfig(
        method=_map_solver_type(str(model.config.get("SolverType", "BDF") or "BDF")),
        rtol=float(rel_tol),
        atol=float(abs_tol),
        max_step=float(max_step),
        seed=seed,
    )
    _log_solver_banner(model, solver_config, days, emit_diagnostics)
    _time_unit_sniff(model, sample_interval_hours)

    if seed is not None:
        np.random.seed(seed)

    state = model.initial_state().astype(float)
    model.apply_initial_assignments_to_state(state)

    if event_log is not None:
        event_log.clear()

    sample_times = _sample_times(days, sample_interval_hours)
    tol_time = _TIME_TOL

    def reconcile(vec: np.ndarray) -> Dict[str, float]:
        ctx = model.build_context_from_state(vec.copy())
        model.evaluate_repeated_assignments(ctx)
        model.apply_algebraic_rules(ctx, vec, mutate=False)
        model.sync_state_from_context(ctx, vec)
        return ctx

    context = reconcile(state)
    _perform_t0_quick_check(model, state, context)

    samples: Dict[float, np.ndarray] = {}
    samples[0.0] = state.copy()

    therapy_flag = therapy.lower()
    base_doses = list(custom_doses) if custom_doses is not None else (
        [] if therapy_flag == "none" else list(model.doses)
    )

    if scheduler is not None:
        scheduled_list = list(scheduler.build(base_doses, model=model, days=days))
    else:
        scheduled_list = _materialise_schedule(base_doses, days)

    if not scheduled_list and therapy_flag == "anti_pd1" and custom_doses is None:
        scheduled_list.extend(_fallback_anti_pd1_doses(model, days))

    scheduled_list.sort(key=lambda item: (item.time, item.priority))
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

    def record_dose_audit(time_days: float, scheduled_dose: ScheduledDose, info: Dict[str, object]) -> None:
        if dose_audit is None:
            return
        dose_audit.append(
            {
                "time_days": time_days,
                "time_hours": time_days * 24.0,
                "dose_name": scheduled_dose.dose.name,
                "target": scheduled_dose.dose.target,
                "interpreted_dimension": info.get("interpreted_dimension"),
                "compartment": info.get("compartment"),
                "compartment_volume_l": info.get("compartment_volume_l"),
                "delta_applied": info.get("delta_applied"),
                "amount_moles": scheduled_dose.amount,
                "amount_mg": scheduled_dose.amount_mg,
                "time_unit": model.time_unit,
            }
        )

    while dose_index < len(scheduled_list) and abs(scheduled_list[dose_index].time) <= tol_time:
        scheduled_dose = scheduled_list[dose_index]
        audit_info = model.apply_dose(scheduled_dose.dose, scheduled_dose.amount, context, current_state)
        context = reconcile(current_state)
        samples[0.0] = current_state.copy()
        record_dose_audit(0.0, scheduled_dose, audit_info)
        logger.info(
            "dose_event time=%g target=%s amount=%g",
            scheduled_dose.time,
            scheduled_dose.dose.target,
            scheduled_dose.amount,
        )
        dose_index += 1

    current_time = 0.0

    while current_time < days - tol_time:
        next_dose_time = scheduled_list[dose_index].time if dose_index < len(scheduled_list) else float("inf")
        next_event_time = pending_events[0].time if pending_events else float("inf")
        target_time = min(days, next_dose_time, next_event_time)

        if target_time <= current_time + tol_time:
            current_time = target_time
            while pending_events and abs(pending_events[0].time - current_time) <= tol_time:
                scheduled_event = heapq.heappop(pending_events)
                for assignment in scheduled_event.entry.assignments:
                    value = assignment.compiled.evaluate(context) if assignment.compiled is not None else 0.0
                    model._apply_target_value(assignment.target, value, context, current_state)
                context = reconcile(current_state)
                samples[current_time] = current_state.copy()
            while dose_index < len(scheduled_list) and scheduled_list[dose_index].time <= current_time + tol_time:
                scheduled_dose = scheduled_list[dose_index]
                audit_info = model.apply_dose(scheduled_dose.dose, scheduled_dose.amount, context, current_state)
                context = reconcile(current_state)
                samples[current_time] = current_state.copy()
                record_dose_audit(current_time, scheduled_dose, audit_info)
                logger.info(
                    "dose_event time=%g target=%s amount=%g",
                    scheduled_dose.time,
                    scheduled_dose.dose.target,
                    scheduled_dose.amount,
                )
                dose_index += 1
            continue

        sol = solve_ivp(
            model.rhs,
            (current_time, target_time),
            current_state,
            method=solver_config.method,
            rtol=solver_config.rtol,
            atol=solver_config.atol,
            max_step=solver_config.max_step,
            dense_output=True,
            events=event_functions if event_functions else None,
            vectorized=False,
        )

        if not sol.success:
            raise NumericsError(f"Integration failed at t={current_time}: {sol.message}")

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
            samples[event_time] = current_state.copy()
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
            samples[current_time] = current_state.copy()

        while dose_index < len(scheduled_list) and scheduled_list[dose_index].time <= current_time + tol_time:
            scheduled_dose = scheduled_list[dose_index]
            audit_info = model.apply_dose(scheduled_dose.dose, scheduled_dose.amount, context, current_state)
            context = reconcile(current_state)
            samples[current_time] = current_state.copy()
            record_dose_audit(current_time, scheduled_dose, audit_info)
            logger.info(
                "dose_event time=%g target=%s amount=%g",
                scheduled_dose.time,
                scheduled_dose.dose.target,
                scheduled_dose.amount,
            )
            dose_index += 1

    if current_time not in samples:
        samples[current_time] = current_state.copy()

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
    contexts: List[Dict[str, float]] = []

    for vector in states:
        vec = vector.copy()
        ctx = model.build_context_from_state(vec)
        model.evaluate_repeated_assignments(ctx)
        model.apply_algebraic_rules(ctx, vec, mutate=False)
        model.sync_state_from_context(ctx, vec)
        contexts.append(dict(ctx))

        c_cells = ctx.get("C1", 0.0)
        d_cells = ctx.get("C_x", 0.0)
        t_tumour = ctx.get("V_T.T1", 0.0)
        tumour_t0 = ctx.get("V_T.T0", 0.0)
        volume_l = ctx.get("V_T", 0.0)
        diameter_cm = ((3.0 * volume_l * 1e3) / (4.0 * math.pi)) ** (1.0 / 3.0) * 2.0 if volume_l > 0 else 0.0
        occupancy = ctx.get("H_PD1_C1", 0.0)
        density = t_tumour / max(volume_l * 1e6, 1e-12)
        t_total = ctx.get("T_total", t_tumour + tumour_t0)

        cancer_cells.append(c_cells)
        dead_cells.append(d_cells)
        t_cells.append(t_total)
        tumour_volume_l.append(volume_l)
        tumour_diameter_cm.append(diameter_cm)
        pd1_occupancy.append(occupancy)
        tcell_density.append(density)

    extras_results: Dict[str, np.ndarray] = {}
    extra_order: List[str] = []

    plugin_sequence: List[ExtraOutputs] = []
    if context_outputs:
        plugin_sequence.append(ContextKeyOutputs(context_outputs))
    if extra_outputs:
        plugin_sequence.extend(extra_outputs)

    if plugin_sequence:
        meta: Dict[str, object] = {
            "solver": solver_config.as_dict(),
            "provenance": model.provenance,
            "time_unit": model.time_unit,
        }
        for plugin in plugin_sequence:
            payload = plugin.compute(sample_times, states, contexts, meta)
            for column, values in payload.items():
                extras_results[column] = np.asarray(values, dtype=float)
                extra_order.append(column)

    header = BASE_HEADER + tuple(extra_order)

    provenance = dict(model.provenance)
    provenance.update(
        {
            "solver_config": json.dumps(solver_config.as_dict(), sort_keys=True),
            "solver_hash": solver_config.identity(),
            "run_label": run_label or "",
        }
    )

    return ScenarioResult(
        time_days=sample_times,
        cancer_cells=np.array(cancer_cells, dtype=float),
        dead_cells=np.array(dead_cells, dtype=float),
        t_cells=np.array(t_cells, dtype=float),
        tumour_volume_l=np.array(tumour_volume_l, dtype=float),
        tumour_diameter_cm=np.array(tumour_diameter_cm, dtype=float),
        pd1_occupancy=np.array(pd1_occupancy, dtype=float),
        tcell_density_per_ul=np.array(tcell_density, dtype=float),
        extras=extras_results,
        header=header,
        provenance=provenance,
        semantics_version=SEMANTICS_VERSION,
    )
