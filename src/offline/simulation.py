"""Numerical runtime for frozen SimBiology snapshots."""

from __future__ import annotations

import heapq
import logging
import math
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .entities import DoseEntry, EventEntry, ScenarioResult, ScheduledDose
from .snapshot import FrozenModel, load_frozen_model, sha256_file, snapshot_digest

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


@dataclass(frozen=True, order=True)
class ScheduledEvent:
    time: float
    priority: int
    entry: EventEntry = field(compare=False)
    trigger_time: float = field(compare=False)
    delay: float = field(compare=False)
    assignments: str = field(compare=False)


# --------------------------------------------------------------------------------------
# Dosing helpers
# --------------------------------------------------------------------------------------


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


# --------------------------------------------------------------------------------------
# Sampling + diagnostics
# --------------------------------------------------------------------------------------


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


# --------------------------------------------------------------------------------------
# Public simulation entry point
# --------------------------------------------------------------------------------------


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
        snapshot_sha = snapshot_digest(model.source_dir)
        equations_sha = sha256_file(equations_path) if equations_path.exists() else "NA"
        config_sha = sha256_file(config_path) if config_path.exists() else "NA"
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

        if not sol.success:  # pragma: no cover
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

    for vector in states:
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
