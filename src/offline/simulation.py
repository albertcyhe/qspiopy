"""Numerical runtime for frozen SimBiology snapshots."""

from __future__ import annotations

import copy
import hashlib
import json
import logging
import math
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Protocol

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from .aliases import inject_output_aliases
from .entities import BASE_HEADER, DoseEntry, EventEntry, ScenarioResult, ScheduledDose, SEMANTICS_VERSION
from .errors import AlignmentFail, ConfigError, NumericsError, SnapshotError
from .initial_conditions import ICOptions, generate_initial_conditions
from .modules import (
    apply_parameter_overrides,
    bind_module_blocks,
    disable_repeated_assignments,
    resolve_module_blocks,
)
from .snapshot import FrozenModel, load_frozen_model
from .units import convert_amount
from .segment_integrator import (
    SolverConfig,
    ScheduledEvent,
    TriggerEventSpec,
    T0Options,
    run_segmented_integration,
    warm_start_quarantine,
)

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
_TIME_KEY_DIGITS = 12
_MAX_SEG_ATTEMPTS = 8
_MICRO_SNAP = 5e-9
_STATE_EVENT_HYSTERESIS = 1e-6
_STATE_EVENT_HYSTERESIS = 1e-6


def _tkey(value: float) -> float:
    """Quantise time points before using as dictionary keys."""
    quantised = round(float(value), _TIME_KEY_DIGITS)
    return float(np.nextafter(quantised, np.inf))


def _merge_close_times(times: np.ndarray) -> np.ndarray:
    """Merge nearly-identical time stamps (≤1 ULP apart)."""
    if times.size == 0:
        return times
    ordered = np.sort(np.asarray(times, dtype=float))
    merged = [ordered[0]]
    for current in ordered[1:]:
        if current <= np.nextafter(merged[-1], np.inf):
            continue
        merged.append(current)
    return np.asarray(merged, dtype=float)


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
    molecular_weight = 1.436e5
    amount_mg = patient_weight * mg_per_kg
    amount = convert_amount(amount_mg, "milligram", molecular_weight)
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
        molecular_weight_g_per_mol=molecular_weight,
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


def _dedupe_scheduled_doses(
    scheduled: Sequence[ScheduledDose],
    tol: float = 1e-12,
) -> List[ScheduledDose]:
    seen: set = set()
    result: List[ScheduledDose] = []
    for item in sorted(scheduled, key=lambda s: (round(s.time, 12), s.priority, s.dose.name)):
        key = (
            round(float(item.time), 12),
            item.dose.name,
            item.dose.target,
            round(float(item.amount or 0.0), 15),
        )
        if key in seen:
            continue
        seen.add(key)
        result.append(item)
    return result


def _guard_single_t0_doses(scheduled: Sequence[ScheduledDose], tol: float = 1e-12) -> List[ScheduledDose]:
    seen: set = set()
    guarded: List[ScheduledDose] = []
    for item in scheduled:
        if abs(item.time) <= tol:
            key = item.dose.target
            if key in seen:
                logger.warning("drop duplicate t=0 dose for target %s", key)
                continue
            seen.add(key)
        guarded.append(item)
    return guarded


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
    candidate_times = _merge_close_times(sample_times)
    for t_sample in candidate_times:
        key = _tkey(t_sample)
        if key in samples:
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
        inject_output_aliases(ctx)
        samples[key] = vec


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


def _effective_max_step(model: FrozenModel, solver: SolverConfig, span: float) -> float:
    if span <= _TIME_TOL:
        return span
    cap = 1.0
    candidate = 0.5 * min(span, cap)
    if candidate <= _TIME_TOL:
        candidate = min(span, cap)
    if math.isfinite(solver.max_step):
        candidate = min(candidate, solver.max_step)
    candidate = min(candidate, span)
    if candidate <= _TIME_TOL:
        candidate = max(span * 0.5, _TIME_TOL * 10.0)
    return candidate


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
    ic_mode: str = "target_volume",
    ic_options: Optional[ICOptions] = None,
    t0_options: Optional[T0Options] = None,
    _disable_warm_start: bool = False,
    param_overrides: Optional[Dict[str, float]] = None,
    module_blocks: Optional[Sequence[str]] = None,
    capture_contexts: bool = False,
) -> ScenarioResult:
    try:
        model = load_frozen_model(snapshot)
    except FileNotFoundError as exc:  # pragma: no cover - passthrough
        raise SnapshotError(str(exc)) from exc
    except Exception as exc:  # pragma: no cover
        raise SnapshotError(f"Failed to load snapshot '{snapshot}': {exc}") from exc

    if param_overrides or module_blocks:
        model = copy.deepcopy(model)

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
    t0_opts = t0_options or T0Options()
    if emit_diagnostics:
        logger.info(
            "t0_options debounce_window=%g bump_eps=%g first_step_cap=%g max_step_cap=%g nan_guard=%s",
            t0_opts.debounce_window_days,
            t0_opts.bump_eps_days,
            t0_opts.first_step_cap_days,
            t0_opts.max_step_cap_days,
            t0_opts.nan_guard,
        )

    if seed is not None:
        np.random.seed(seed)

    if param_overrides:
        apply_parameter_overrides(model, dict(param_overrides))

    if module_blocks is None:
        block_names: Tuple[str, ...] = ("tumour_geometry", "pd1_bridge")
    else:
        block_names = tuple(module_blocks)
    pre_blocks, post_blocks, disable_from_blocks, active_block_names = resolve_module_blocks(
        model, block_names
    )
    if pre_blocks or post_blocks:
        bind_module_blocks(model, pre_blocks, post_blocks)
    disable_targets = tuple(dict.fromkeys(disable_from_blocks))
    if disable_targets:
        disable_repeated_assignments(model, disable_targets)
    if active_block_names:
        logger.info("module_blocks active=%s", ",".join(active_block_names))

    ic_fallback = False
    if ic_mode == "target_volume":
        ic_opts = ic_options or ICOptions(
            target_diameter_cm=0.5,
            max_days=150.0,
            max_wall_seconds=20.0,
            reset_policy="cancer_only",
        )
        try:
            state, context = generate_initial_conditions(model, opts=ic_opts)
            inject_output_aliases(context)
        except ValueError as exc:
            logger.warning("target_volume IC not available for snapshot '%s' (%s); falling back to snapshot IC", snapshot, exc)
            ic_fallback = True
    if ic_mode != "target_volume" or ic_fallback:
        state = model.initial_state().astype(float)
        model.apply_initial_assignments_to_state(state)
        context = None

    if event_log is not None:
        event_log.clear()

    sample_times = _sample_times(days, sample_interval_hours)
    sample_times = _merge_close_times(sample_times)
    tol_time = _TIME_TOL

    def _reconcile_mutating(vec: np.ndarray, time_point: Optional[float] = None) -> Dict[str, float]:
        ctx = model.build_context_from_state(vec.copy())
        if time_point is not None:
            t_val = float(time_point)
            ctx["time"] = t_val
            ctx["time_days"] = t_val
            ctx["t"] = t_val
        model.evaluate_repeated_assignments(ctx)
        model.apply_algebraic_rules(ctx, vec, mutate=False)
        model.sync_state_from_context(ctx, vec)
        inject_output_aliases(ctx)
        return ctx

    def _reconcile_readonly(vec: np.ndarray, time_point: Optional[float] = None) -> Dict[str, float]:
        ctx = model.build_context_from_state(vec.copy())
        if time_point is not None:
            t_val = float(time_point)
            ctx["time"] = t_val
            ctx["time_days"] = t_val
            ctx["t"] = t_val
        model.evaluate_repeated_assignments(ctx)
        model.apply_algebraic_rules(ctx, vec, mutate=False)
        inject_output_aliases(ctx)
        return ctx

    context = _reconcile_mutating(state, time_point=0.0)
    _perform_t0_quick_check(model, state, context)

    samples: Dict[float, np.ndarray] = {}
    samples[_tkey(0.0)] = state.copy()

    therapy_flag = therapy.lower()
    base_doses = list(custom_doses) if custom_doses is not None else (
        [] if therapy_flag == "none" else list(model.doses)
    )

    if scheduler is not None:
        scheduled_list = list(scheduler.build(base_doses, model=model, days=days))
    else:
        scheduled_list = _materialise_schedule(base_doses, days)
    scheduled_list = _dedupe_scheduled_doses(scheduled_list)
    scheduled_list = _guard_single_t0_doses(scheduled_list)

    if not scheduled_list and therapy_flag == "anti_pd1" and custom_doses is None:
        existing_signatures = {
            (dose.dose.name, dose.dose.target, round(dose.time, 6), round(dose.dose.amount or 0.0, 12))
            for dose in scheduled_list
        }
        fallback = _fallback_anti_pd1_doses(model, days)
        filtered = [
            dose
            for dose in fallback
            if (dose.dose.name, dose.dose.target, round(dose.time, 6), round(dose.dose.amount or 0.0, 12))
            not in existing_signatures
        ]
        scheduled_list.extend(filtered)

    scheduled_list.sort(key=lambda item: (item.time, item.priority))
    if emit_diagnostics:
        logger.info("scheduled_doses=%d", len(scheduled_list))
    model._active_dose_schedule = tuple(scheduled_list)
    dose_index = 0
    pending_events: List[ScheduledEvent] = []
    current_state = state

    def make_trigger_spec(entry: EventEntry) -> TriggerEventSpec:
        def fn(t: float, y: np.ndarray) -> float:
            vec = np.asarray(y, dtype=float)
            ctx = model.build_context_from_state(vec.copy())
            model.evaluate_repeated_assignments(ctx)
            model.apply_algebraic_rules(ctx, vec, mutate=False)
            return float(entry.trigger_compiled.evaluate(ctx))

        fn.direction = entry.direction
        fn.terminal = False
        return TriggerEventSpec(entry=entry, fn=fn, state={})

    trigger_specs = [make_trigger_spec(entry) for entry in model.events]

    phase_codes = {"cont": 0, "pre": 1, "post": 2}
    kind_codes = {"cont": 0, "dose": 1, "event": 2}
    record_annotations: Dict[str, Dict[float, object]] = {
        "phase_code": {},
        "discontinuity_type": {},
        "event_index": {},
        "event_name": {},
        "dose_name": {},
        "dose_target": {},
    }

    def record_state_annotation(
        time_days: float,
        vector: np.ndarray,
        phase: str,
        kind: str,
        event_index: Optional[int],
        event_name: Optional[str],
        dose_target: Optional[str],
        dose_name: Optional[str],
    ) -> None:
        key = _tkey(time_days)
        samples[key] = vector.copy()
        record_annotations["phase_code"][key] = phase_codes.get(phase, 0)
        record_annotations["discontinuity_type"][key] = kind_codes.get(kind, 0)
        if event_index is not None:
            record_annotations["event_index"][key] = event_index
        if event_name:
            record_annotations["event_name"][key] = event_name
        if dose_name:
            record_annotations["dose_name"][key] = dose_name
        if dose_target:
            record_annotations["dose_target"][key] = dose_target

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
                "units": info.get("units"),
                "compartment": info.get("compartment"),
                "compartment_volume_l": info.get("compartment_volume_l"),
                "delta_state_value": info.get("delta_state_value"),
                "delta_amount_mol": info.get("delta_amount_mol"),
                "amount_moles": scheduled_dose.amount,
                "amount_mg": scheduled_dose.amount_mg,
                "time_unit": model.time_unit,
            }
        )

    while dose_index < len(scheduled_list) and abs(scheduled_list[dose_index].time) <= tol_time:
        scheduled_dose = scheduled_list[dose_index]
        audit_info = model.apply_dose(scheduled_dose.dose, scheduled_dose.amount, context, current_state)
        context = _reconcile_mutating(current_state, time_point=0.0)
        samples[_tkey(0.0)] = current_state.copy()
        record_dose_audit(0.0, scheduled_dose, audit_info)
        logger.info(
            "dose_event time=%g target=%s amount=%g",
            scheduled_dose.time,
            scheduled_dose.dose.target,
            scheduled_dose.amount,
        )
        dose_index += 1

    jac_pattern = model.jacobian_sparsity()

    current_state = state
    if not _disable_warm_start:
        warm_state, context, start_time = warm_start_quarantine(
            model=model,
            solver_config=solver_config,
            state=current_state,
            context=context,
            reconcile=_reconcile_mutating,
            samples=samples,
            time_key=_tkey,
            pending_events=pending_events,
            t0_options=t0_opts,
            jac_sparsity=jac_pattern,
            tol_time=tol_time,
            diagnostics=emit_diagnostics,
        )
        current_state = warm_state
        samples.setdefault(_tkey(start_time), current_state.copy())
    else:
        start_time = 0.0
        samples.setdefault(_tkey(start_time), current_state.copy())
        context = _reconcile_mutating(current_state, time_point=start_time)
    initial_trigger_bools = [_evaluate_trigger_bool(entry, context) for entry in model.events]

    current_state, context, current_time, dose_index = run_segmented_integration(
        model=model,
        solver_config=solver_config,
        scheduled_doses=scheduled_list,
        pending_events=pending_events,
        trigger_specs=trigger_specs,
        initial_trigger_bools=initial_trigger_bools,
        sample_times=sample_times,
        samples=samples,
        state=current_state,
        context=context,
        start_time=start_time,
        stop_time=days,
        dose_index=dose_index,
        tol_time=tol_time,
        record_solution_samples=_record_solution_samples,
        record_dose_audit=record_dose_audit,
        reconcile=_reconcile_mutating,
        time_key=_tkey,
        event_log=event_log,
        jac_sparsity=jac_pattern,
        record_state=record_state_annotation,
        t0_options=t0_opts,
        diagnostics=emit_diagnostics,
    )

    sample_key_set = {_tkey(value) for value in sample_times}
    recorded_keys = sample_key_set.intersection(samples.keys())
    expected_count = len(sample_key_set)
    recorded_count = len(recorded_keys)
    if expected_count:
        min_required = min(expected_count, max(5, int(0.8 * expected_count)))
    else:
        min_required = 0
    if emit_diagnostics and expected_count:
        missing_keys = sorted(sample_key_set.difference(samples.keys()))
        sample_preview = sorted(samples.keys())[:5]
        logger.warning(
            "sample_capture expected=%d recorded=%d missing_first=%s have_first=%s",
            expected_count,
            recorded_count,
            missing_keys[:5],
            sample_preview,
        )

    if not _disable_warm_start and expected_count and recorded_count < min_required:
        logger.warning(
            "warm-start captured %d/%d sample states; retrying without warm-start",
            recorded_count,
            expected_count,
        )
        return simulate_frozen_model(
            snapshot,
            days=days,
            therapy=therapy,
            rtol_override=rtol_override,
            atol_override=atol_override,
            seed=seed,
            sample_interval_hours=sample_interval_hours,
            run_label=run_label,
            scheduler=scheduler,
            custom_doses=custom_doses,
            extra_outputs=extra_outputs,
            event_log=event_log,
            dose_audit=dose_audit,
            context_outputs=context_outputs,
            ic_mode=ic_mode,
            ic_options=ic_options,
            t0_options=t0_options,
            param_overrides=param_overrides,
            module_blocks=module_blocks,
            emit_diagnostics=emit_diagnostics,
            _disable_warm_start=True,
        )

    current_key = _tkey(current_time)
    if current_key not in samples:
        samples[current_key] = current_state.copy()

    output_states = []
    last_state = current_state.copy()
    for t_val in sample_times:
        key = _tkey(t_val)
        if key in samples:
            last_state = samples[key]
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
    tcell_params = model.parameters or {}
    tcell_tau_days = max(float(tcell_params.get("tcell_alignment_tau_days", 0.0)), 0.0)
    tcell_w_live = float(tcell_params.get("tcell_alignment_w_live", 1.0))
    tcell_w_treg = float(tcell_params.get("tcell_alignment_w_treg", 1.0))
    tcell_offset = float(tcell_params.get("tcell_alignment_offset_cells", 0.0))
    tcell_min_cells = max(float(tcell_params.get("tcell_alignment_min_cells", 0.0)), 0.0)
    tcell_occ_supp = float(tcell_params.get("tcell_alignment_occ_supp_coeff", 0.0))
    tcell_state_output: Optional[float] = None
    prev_time_val: Optional[float] = None

    for time_val, vector in zip(sample_times, states):
        vec = vector.copy()
        ctx = _reconcile_readonly(vec, time_point=time_val)
        contexts.append(dict(ctx))

        c_cells = float(ctx.get("C1", 0.0))
        d_cells = float(ctx.get("C_x", 0.0))
        t_tumour = float(ctx.get("V_T.T1", 0.0))
        tumour_t0 = float(ctx.get("V_T.T0", 0.0))
        volume_l = float(ctx.get("tumour_volume_l", ctx.get("V_T", 0.0)))
        diameter_cm = float(
            ctx.get(
                "tumour_diameter_cm",
                ((3.0 * volume_l * 1e3) / (4.0 * math.pi)) ** (1.0 / 3.0) * 2.0 if volume_l > 0.0 else 0.0,
            )
        )
        occupancy = float(ctx.get("pd1_occupancy", ctx.get("H_PD1_C1", 0.0)))
        density = float(ctx.get("tcell_density_per_ul", 0.0))
        if not density and volume_l > 0.0:
            density = t_tumour / max(volume_l * 1e6, 1e-12)

        target_tcells = max(
            tcell_offset + tcell_w_live * t_tumour + tcell_w_treg * tumour_t0,
            tcell_min_cells,
        )
        if tcell_occ_supp:
            target_tcells *= max(0.0, 1.0 - tcell_occ_supp * occupancy)
        delta_t = 0.0 if prev_time_val is None else max(0.0, float(time_val - prev_time_val))
        prev_time_val = float(time_val)

        if tcell_state_output is None:
            tcell_state_output = target_tcells
        elif tcell_tau_days > 0.0 and delta_t > 0.0:
            beta = math.exp(-delta_t / tcell_tau_days)
            tcell_state_output = target_tcells + (tcell_state_output - target_tcells) * beta
        else:
            tcell_state_output = target_tcells
        t_total = tcell_state_output

        cancer_cells.append(c_cells)
        dead_cells.append(d_cells)
        t_cells.append(t_total)
        tumour_volume_l.append(volume_l)
        tumour_diameter_cm.append(diameter_cm)
        pd1_occupancy.append(occupancy)
        tcell_density.append(density)

    raw_contexts_tuple: Tuple[Dict[str, float], ...] = tuple(contexts) if capture_contexts else tuple()
    raw_states = states.copy() if capture_contexts else None

    extras_results: Dict[str, np.ndarray] = {}
    base_extras = [
        ("phase_code", record_annotations["phase_code"], float, 0.0),
        ("discontinuity_type_code", record_annotations["discontinuity_type"], float, 0.0),
        ("event_index", record_annotations["event_index"], float, 0.0),
        ("event_name", record_annotations["event_name"], object, ""),
        ("dose_name", record_annotations["dose_name"], object, ""),
        ("dose_target", record_annotations["dose_target"], object, ""),
    ]
    extra_order: List[str] = []
    for column, mapping, dtype, default in base_extras:
        extras_results[column] = np.array(
            [mapping.get(_tkey(t_val), default) for t_val in sample_times],
            dtype=dtype,
        )
        extra_order.append(column)

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
        raw_states=raw_states,
        raw_contexts=raw_contexts_tuple,
    )
def _evaluate_trigger_bool(entry: EventEntry, context: Dict[str, float]) -> bool:
    if entry.trigger_boolean_compiled is not None:
        return bool(entry.trigger_boolean_compiled.evaluate(context))
    value = entry.trigger_compiled.evaluate(context)
    return float(value) >= 0.0
