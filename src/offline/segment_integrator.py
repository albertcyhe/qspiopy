"""Segmented integration helper for SimBiology-derived models."""

from __future__ import annotations

import hashlib
import heapq
import json
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .entities import EventEntry, ScheduledDose
from .errors import NumericsError
from .snapshot import FrozenModel
import time

_MICRO_SNAP = 5e-9
_STATE_EVENT_SUPPRESS_WINDOW = 1e-6
_EVENT_HYSTERESIS = 1e-6


@dataclass
class DebounceOptions:
    refractory_dt: float = 1e-9
    rearm_requires_false: bool = True
    rearm_eps: float = 1e-9
    debounce_window_days: float = 1e-8


@dataclass(frozen=True)
class T0Options:
    """Options that govern warm-start behaviour and early-time damping."""

    rearm_eps: float = 1e-9
    debounce_window_days: float = 1e-8
    bump_eps_days: float = 1e-8
    first_step_cap_days: float = 1e-6
    max_step_cap_days: float = 1e-4
    tiny_segment_trials: int = 6
    tiny_segment_wall_s: Optional[float] = 2.0
    nan_guard: bool = True
    refractory_dt: float = 1e-9


class EventDebouncer:
    def __init__(self, initial_bools: Sequence[bool], opts: DebounceOptions, start_time: float = 0.0):
        self.opts = opts
        self.armed = [not bool(val) for val in initial_bools]
        self.refractory_until = [-float("inf")] * len(initial_bools)
        window = max(float(opts.debounce_window_days or 0.0), 0.0)
        self.window_end = start_time + window

    def is_armed(self, idx: int, t: float) -> bool:
        return self.armed[idx] and t >= self.refractory_until[idx]

    def should_track(self, idx: int, t: float) -> bool:
        if float(t) < self.window_end:
            return False
        return self.is_armed(idx, float(t))

    def fire(self, idx: int, t: float) -> None:
        self.armed[idx] = False
        self.refractory_until[idx] = t + self.opts.refractory_dt

    def observe(self, idx: int, t: float, trigger_value: float) -> None:
        is_true = trigger_value >= 0.0
        if self.armed[idx]:
            return
        if t < self.refractory_until[idx]:
            return
        if not self.opts.rearm_requires_false:
            self.armed[idx] = True
            return
        if not is_true:
            threshold = -abs(self.opts.rearm_eps or 0.0)
            if trigger_value <= threshold:
                self.armed[idx] = True


@dataclass(frozen=True)
class SolverConfig:
    """Configuration driving scipy's solve_ivp."""

    method: str
    rtol: float
    atol: float
    max_step: float
    seed: Optional[int] = None

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


@dataclass(frozen=True, order=True)
class ScheduledEvent:
    time: float
    priority: int
    entry: EventEntry = field(compare=False)
    trigger_time: float = field(compare=False)
    delay: float = field(compare=False)
    assignments: str = field(compare=False)


@dataclass(frozen=True)
class TriggerEventSpec:
    entry: EventEntry
    fn: Callable[[float, np.ndarray], float]
    state: Dict[str, bool]


def _process_pending_events(
    *,
    time_point: float,
    pending_events: List[ScheduledEvent],
    model: FrozenModel,
    context: Dict[str, float],
    state: np.ndarray,
    reconcile: Callable[[np.ndarray, float], Dict[str, float]],
    debouncer: EventDebouncer,
    event_positions: Dict[int, int],
    samples: Dict[float, np.ndarray],
    time_key: Callable[[float], float],
    event_log: Optional[List[Dict[str, object]]],
    tol_time: float,
) -> bool:
    fired = False
    while pending_events and abs(pending_events[0].time - time_point) <= tol_time:
        fired = True
        scheduled_event = heapq.heappop(pending_events)
        if event_log is not None:
            event_log.append(
                {
                    "event_index": scheduled_event.entry.index,
                    "time_fire": time_point,
                    "time_trigger": scheduled_event.trigger_time,
                    "delay": scheduled_event.delay,
                    "type": "delayed",
                    "assignments": scheduled_event.assignments,
                }
            )
        for assignment in scheduled_event.entry.assignments:
            value = assignment.compiled.evaluate(context) if assignment.compiled is not None else 0.0
            model._apply_target_value(assignment.target, value, context, state)
        pos = event_positions.get(scheduled_event.entry.index)
        if pos is not None:
            debouncer.fire(pos, time_point)
        context = reconcile(state, time_point)
        samples[time_key(time_point)] = state.copy()
    return fired


def _coalesce_same_time_doses(
    doses: Sequence[ScheduledDose],
    tol: float = 1e-12,
) -> List[ScheduledDose]:
    buckets: Dict[Tuple[str, str], Dict[str, float]] = {}
    exemplars: Dict[Tuple[str, str], ScheduledDose] = {}
    for dose in doses:
        key = (dose.dose.name, dose.dose.target)
        info = buckets.setdefault(key, {"amount": 0.0, "amount_mg": 0.0})
        info["amount"] += float(dose.amount or 0.0)
        if dose.amount_mg is not None:
            info["amount_mg"] += float(dose.amount_mg)
        exemplars.setdefault(key, dose)
    collapsed: List[ScheduledDose] = []
    for key, totals in buckets.items():
        exemplar = exemplars[key]
        collapsed.append(
            ScheduledDose(
                time=exemplar.time,
                priority=exemplar.priority,
                dose=exemplar.dose,
                amount=totals["amount"],
                amount_mg=totals["amount_mg"] if totals["amount_mg"] else None,
            )
        )
    collapsed.sort(key=lambda d: (d.priority, d.dose.name))
    return collapsed


def _trigger_value(entry: EventEntry, context: Dict[str, float]) -> float:
    if entry.trigger_compiled is not None:
        return float(entry.trigger_compiled.evaluate(context))
    if entry.trigger_boolean_compiled is not None:
        raw = entry.trigger_boolean_compiled.evaluate(context)
        return 1.0 if bool(raw) else -1.0
    return -1.0


def _trigger_bool(entry: EventEntry, context: Dict[str, float]) -> bool:
    if entry.trigger_boolean_compiled is not None:
        return bool(entry.trigger_boolean_compiled.evaluate(context))
    value = _trigger_value(entry, context)
    return value >= 0.0


def _bounded_first_step(current: float, limit: float, cap: Optional[float], enforce: bool) -> float:
    value = _initial_step(current, limit)
    if not enforce:
        return value
    if cap is None or cap <= 0.0:
        return value
    return min(value, cap)


def _bounded_max_step(solver_cap: float, cap_override: Optional[float], enforce: bool) -> float:
    base = solver_cap
    if not enforce or cap_override is None or cap_override <= 0.0:
        return base
    if np.isfinite(base):
        return min(base, cap_override)
    return cap_override


def _nan_guard(model: FrozenModel, state: np.ndarray, t_point: float) -> None:
    values = model.rhs(t_point, state.copy())
    if not np.all(np.isfinite(values)):
        raise NumericsError("t=0 warm-start: non-finite derivative detected")


def _warm_start_segment(
    *,
    model: FrozenModel,
    solver_config: SolverConfig,
    state: np.ndarray,
    context: Dict[str, float],
    start_time: float,
    reconcile: Callable[[np.ndarray, float], Dict[str, float]],
    samples: Dict[float, np.ndarray],
    time_key: Callable[[float], float],
    t0_options: T0Options,
    tol_time: float,
    jac_sparsity: Optional[np.ndarray],
) -> Tuple[np.ndarray, Dict[str, float], float]:
    if t0_options.nan_guard:
        _nan_guard(model, state, start_time)

    span = max(float(t0_options.bump_eps_days or 0.0), 0.0)
    if span <= tol_time:
        return state, context, start_time

    attempts = max(int(t0_options.tiny_segment_trials or 1), 1)
    wall_limit = t0_options.tiny_segment_wall_s
    start_wall = time.perf_counter()
    target_time = start_time + span
    for attempt in range(attempts):
        first_step = min(span, t0_options.first_step_cap_days or span)
        if first_step <= 0.0:
            first_step = span
        max_step = _bounded_max_step(solver_config.max_step, t0_options.max_step_cap_days, True)
        max_step = min(max_step, span) if np.isfinite(max_step) else span
        sol = solve_ivp(
            model.rhs,
            (start_time, target_time),
            state.copy(),
            method=solver_config.method,
            rtol=solver_config.rtol,
            atol=solver_config.atol,
            first_step=first_step,
            max_step=max_step,
            dense_output=False,
            events=None,
            jac_sparsity=jac_sparsity,
            vectorized=False,
        )
        if sol.success and sol.y.size:
            warm_state = np.asarray(sol.y[:, -1], dtype=float)
            context = reconcile(warm_state, target_time)
            samples.setdefault(time_key(target_time), warm_state.copy())
            return warm_state, context, target_time
        if wall_limit is not None and (time.perf_counter() - start_wall) >= wall_limit:
            break
        span *= 0.5
        if span <= tol_time:
            break
        target_time = start_time + span
    raise NumericsError("Warm-start segment failed to advance away from t=0")


def run_segmented_integration(
    *,
    model: FrozenModel,
    solver_config: SolverConfig,
    scheduled_doses: Sequence[ScheduledDose],
    pending_events: List[ScheduledEvent],
    trigger_specs: Sequence[TriggerEventSpec],
    initial_trigger_bools: Sequence[bool],
    sample_times: np.ndarray,
    samples: Dict[float, np.ndarray],
    state: np.ndarray,
    context: Dict[str, float],
    start_time: float,
    stop_time: float,
    dose_index: int,
    tol_time: float,
    record_solution_samples: Callable[..., None],
    record_dose_audit: Callable[[float, ScheduledDose, Dict[str, object]], None],
    reconcile: Callable[[np.ndarray], Dict[str, float]],
    time_key: Callable[[float], float],
    event_log: Optional[List[Dict[str, object]]] = None,
    jac_sparsity: Optional[np.ndarray] = None,
    record_state: Optional[
        Callable[
            [float, np.ndarray, str, str, Optional[int], Optional[str], Optional[str], Optional[str]],
            None,
        ]
    ] = None,
    t0_options: Optional[T0Options] = None,
) -> Tuple[np.ndarray, Dict[str, float], float, int]:
    """Integrate between discontinuities while replaying events/doses."""

    segment_points = sorted(
        {
            start_time,
            stop_time,
            *[
                min(stop_time, scheduled_dose.time)
                for scheduled_dose in scheduled_doses
                if scheduled_dose.time <= stop_time + tol_time
            ],
        }
    )
    options = t0_options or T0Options()
    current_time = start_time
    segment_index = 0
    event_positions = {spec.entry.index: idx for idx, spec in enumerate(trigger_specs)}
    debounce_config = DebounceOptions(
        refractory_dt=options.refractory_dt,
        rearm_requires_false=True,
        rearm_eps=options.rearm_eps,
        debounce_window_days=options.debounce_window_days,
    )
    debouncer = EventDebouncer(initial_trigger_bools, debounce_config, start_time=start_time)
    limit_first_step = True
    limit_max_step = True

    trigger_functions: List[Callable[[float, np.ndarray], float]] = []
    for spec in trigger_specs:
        pos = event_positions[spec.entry.index]
        original_fn = spec.fn

        def wrapped(t: float, y: np.ndarray, _orig=original_fn, _pos=pos):
            if not debouncer.should_track(_pos, float(t)):
                return -1.0
            return _orig(t, y)

        wrapped.direction = getattr(original_fn, "direction", spec.entry.direction)
        wrapped.terminal = getattr(original_fn, "terminal", False)
        trigger_functions.append(wrapped)

    trigger_entries = [spec.entry for spec in trigger_specs]
    suppressed_until: Dict[int, float] = {}

    def refresh_debouncer(time_point: float, ctx: Dict[str, float]) -> None:
        for spec in trigger_specs:
            pos = event_positions.get(spec.entry.index)
            if pos is None:
                continue
            trigger_val = _trigger_value(spec.entry, ctx)
            debouncer.observe(pos, time_point, trigger_val)

    def reconcile_with_refresh(vec: np.ndarray, t_point: float) -> Dict[str, float]:
        ctx_local = reconcile(vec)
        refresh_debouncer(t_point, ctx_local)
        return ctx_local

    context = reconcile_with_refresh(state, start_time)
    state, context, current_time = _warm_start_segment(
        model=model,
        solver_config=solver_config,
        state=state,
        context=context,
        start_time=start_time,
        reconcile=reconcile_with_refresh,
        samples=samples,
        time_key=time_key,
        t0_options=options,
        tol_time=tol_time,
        jac_sparsity=jac_sparsity,
    )

    while segment_index < len(segment_points) - 1 and current_time < stop_time - tol_time:
        target_time = segment_points[segment_index + 1]
        reached_target = False

        while not reached_target:
            guard_window = max(1e-4, 10.0 * tol_time)
            remaining = target_time - current_time
            if remaining <= guard_window:
                samples.setdefault(time_key(target_time), state.copy())
                current_time = float(np.nextafter(target_time, np.inf))
                reached_target = True
                break
            segment_end = target_time
            if pending_events:
                segment_end = min(segment_end, pending_events[0].time)
            integration_stop = min(segment_end, target_time - guard_window)
            if integration_stop <= current_time + tol_time:
                reached_target = abs(segment_end - target_time) <= tol_time
                break

            first_step = _bounded_first_step(
                current_time,
                integration_stop,
                options.first_step_cap_days,
                limit_first_step,
            )
            max_step = _bounded_max_step(
                solver_config.max_step,
                options.max_step_cap_days,
                limit_max_step,
            )
            sol = solve_ivp(
                model.rhs,
                (current_time, integration_stop),
                state,
                method=solver_config.method,
                rtol=solver_config.rtol,
                atol=solver_config.atol,
                max_step=max_step,
                first_step=first_step,
                dense_output=True,
                events=trigger_functions if trigger_functions else None,
                jac_sparsity=jac_sparsity,
                vectorized=False,
            )

            if not sol.success:
                message = (sol.message or "").lower()
                if "step size" in message or "strictly increasing" in message:
                    samples.setdefault(time_key(integration_stop), state.copy())
                    current_time = float(np.nextafter(integration_stop, np.inf))
                    break
                raise NumericsError(f"Integration failed at t={current_time}: {sol.message}")

            limit_first_step = False
            limit_max_step = False

            triggered_pairs: List[Tuple[float, EventEntry]] = []
            if sol.t_events and trigger_functions:
                for idx, times in enumerate(sol.t_events):
                    entry = trigger_entries[idx]
                    for time_val in np.atleast_1d(times):
                        triggered_pairs.append((float(time_val), entry))

            if triggered_pairs:
                triggered_pairs.sort(key=lambda pair: (pair[0], pair[1].index))
                event_time = triggered_pairs[0][0]
                if record_state:
                    record_state(event_time, state.copy(), "pre", "event", triggered_pairs[0][1].index, triggered_pairs[0][1].name, None, None)
                record_solution_samples(samples, sample_times, current_time, event_time, sol, model, inclusive=False)
                if sol.sol is not None:
                    state = np.asarray(sol.sol(event_time), dtype=float)
                else:
                    idx = int(np.argmin(np.abs(sol.t - event_time)))
                    state = np.asarray(sol.y[:, idx], dtype=float)
                context = reconcile_with_refresh(state, event_time)
                same_time_entries = [
                    entry for time_val, entry in triggered_pairs if abs(time_val - event_time) <= tol_time
                ]
                same_time_entries = [
                    entry
                    for entry in same_time_entries
                    if suppressed_until.get(entry.index, -np.inf) <= event_time - tol_time
                ]
                if not same_time_entries:
                    current_time = float(np.nextafter(event_time, np.inf))
                    continue
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
                    pos = event_positions.get(entry.index)
                    if pos is not None:
                        debouncer.fire(pos, event_time)
                    for assignment in entry.assignments:
                        value = assignment.compiled.evaluate(context) if assignment.compiled is not None else 0.0
                        model._apply_target_value(assignment.target, value, context, state)
                    context = reconcile_with_refresh(state, event_time)
                    suppressed_until[entry.index] = event_time + max(_STATE_EVENT_SUPPRESS_WINDOW, 10.0 * tol_time)
                samples[time_key(event_time)] = state.copy()
                if record_state:
                    last_entry = same_time_entries[-1] if same_time_entries else triggered_pairs[0][1]
                    record_state(event_time, state.copy(), "post", "event", last_entry.index, last_entry.name, None, None)
                limit_first_step = True
                limit_max_step = True
                current_time = float(np.nextafter(event_time, np.inf))
                continue

            record_solution_samples(
                samples,
                sample_times,
                current_time,
                integration_stop,
                sol,
                model,
                inclusive=True,
            )
            state = np.asarray(sol.y[:, -1], dtype=float)
            context = reconcile_with_refresh(state, current_time)
            current_time = integration_stop

            pending_triggered = _process_pending_events(
                time_point=current_time,
                pending_events=pending_events,
                model=model,
                context=context,
                state=state,
                reconcile=reconcile_with_refresh,
                debouncer=debouncer,
                event_positions=event_positions,
                samples=samples,
                time_key=time_key,
                event_log=event_log,
                tol_time=tol_time,
            )
            if pending_triggered:
                limit_first_step = True
                limit_max_step = True

            if abs(current_time - target_time) <= tol_time:
                reached_target = True

        current_time = target_time

        pending_triggered = _process_pending_events(
            time_point=current_time,
            pending_events=pending_events,
            model=model,
            context=context,
            state=state,
            reconcile=reconcile_with_refresh,
            debouncer=debouncer,
            event_positions=event_positions,
            samples=samples,
            time_key=time_key,
            event_log=event_log,
            tol_time=tol_time,
        )
        if pending_triggered:
            limit_first_step = True
            limit_max_step = True

        dose_applied = False
        while dose_index < len(scheduled_doses) and abs(scheduled_doses[dose_index].time - current_time) <= tol_time:
            batch: List[ScheduledDose] = []
            target_time = scheduled_doses[dose_index].time
            while dose_index < len(scheduled_doses) and abs(scheduled_doses[dose_index].time - target_time) <= tol_time:
                batch.append(scheduled_doses[dose_index])
                dose_index += 1
            for scheduled_dose in _coalesce_same_time_doses(batch, tol=tol_time):
                if record_state:
                    record_state(current_time, state.copy(), "pre", "dose", None, scheduled_dose.dose.name, scheduled_dose.dose.target, scheduled_dose.dose.name)
                audit_info = model.apply_dose(scheduled_dose.dose, scheduled_dose.amount, context, state)
                context = reconcile_with_refresh(state, current_time)
                samples[time_key(current_time)] = state.copy()
                record_dose_audit(current_time, scheduled_dose, audit_info)
                if record_state:
                    record_state(current_time, state.copy(), "post", "dose", None, scheduled_dose.dose.name, scheduled_dose.dose.target, scheduled_dose.dose.name)
                dose_applied = True
        if dose_applied:
            limit_first_step = True
            limit_max_step = True
        if current_time >= stop_time - tol_time:
            break

        current_time = float(np.nextafter(current_time, np.inf))
        segment_index += 1

    return state, context, current_time, dose_index
def _initial_step(current: float, limit: float) -> float:
    span = max(limit - current, 0.0)
    if span <= 0.0:
        return np.finfo(float).eps
    candidate = min(span * 0.05, 0.1)
    floor = np.finfo(float).eps * max(1.0, abs(current))
    return max(candidate, floor)
