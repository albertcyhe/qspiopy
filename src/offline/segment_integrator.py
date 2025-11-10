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
import logging

logger = logging.getLogger(__name__)

_MICRO_SNAP = 5e-9
_STATE_EVENT_SUPPRESS_WINDOW = 1e-6
_EVENT_HYSTERESIS = 1e-6


@dataclass
class DebounceOptions:
    refractory_dt: float = 1e-9
    rearm_requires_false: bool = True
    rearm_eps: float = 1e-9
    rearm_below: Optional[float] = None
    debounce_window_days: float = 1e-8
    freeze_window_days: float = 1e-8
    min_quiet_days: float = 1e-8


@dataclass(frozen=True)
class T0Options:
    """Options that govern warm-start behaviour and early-time damping."""

    rearm_eps: float = 1e-9
    debounce_window_days: float = 1e-8
    freeze_window_days: float = 1e-8
    settle_no_event_until_days: float = 0.0
    bump_eps_days: float = 1e-8
    first_step_cap_days: float = 1e-6
    max_step_cap_days: float = 1e-4
    nan_guard: bool = True
    refractory_dt: float = 1e-9
    quarantine_days: float = 1e-7
    tiny_segment_trials: int = 40
    max_trials: int = 20
    max_wall_seconds: float = 2.0
    rearm_below: Optional[float] = None
    disable_state_events_during_warm: bool = True
    fire_t0_true_events_once: bool = True


class EventDebouncer:
    def __init__(self, initial_bools: Sequence[bool], opts: DebounceOptions, start_time: float = 0.0):
        self.opts = opts
        self.armed = [not bool(val) for val in initial_bools]
        self.refractory_until = [-float("inf")] * len(initial_bools)
        self.last_false_since = [float("inf")] * len(initial_bools)
        window = max(float(opts.debounce_window_days or 0.0), 0.0)
        self.window_end = start_time + window
        freeze = max(float(opts.freeze_window_days or 0.0), 0.0)
        self.freeze_until = start_time + freeze
        self.min_quiet = max(float(opts.min_quiet_days or 0.0), 0.0)
        threshold = opts.rearm_below
        if threshold is None:
            threshold = -abs(opts.rearm_eps or 0.0)
        self.rearm_threshold = float(threshold)

    def fire(self, idx: int, t: float) -> None:
        self.armed[idx] = False
        self.refractory_until[idx] = t + self.opts.refractory_dt
        self.last_false_since[idx] = float("inf")

    def preset_disarm(self, idx: int) -> None:
        if 0 <= idx < len(self.armed):
            self.armed[idx] = False
            self.refractory_until[idx] = -float("inf")
            self.last_false_since[idx] = float("inf")

    def observe(self, idx: int, t: float, trigger_value: float, allow_fire: bool = True) -> float:
        """Return value fed to solver while updating hysteresis state."""
        threshold = self.rearm_threshold
        if float(t) < self.freeze_until:
            allow_fire = False
        ready = self.armed[idx] and t >= self.refractory_until[idx]
        if trigger_value < threshold:
            if self.last_false_since[idx] == float("inf"):
                self.last_false_since[idx] = float(t)
        else:
            self.last_false_since[idx] = float("inf")

        if not ready and t >= self.refractory_until[idx]:
            if self.last_false_since[idx] != float("inf"):
                quiet_elapsed = float(t) - self.last_false_since[idx]
            else:
                quiet_elapsed = 0.0
            if ((not self.opts.rearm_requires_false) or (trigger_value <= threshold)) and (
                quiet_elapsed >= self.min_quiet
            ):
                self.armed[idx] = True
                self.last_false_since[idx] = float("inf")
                ready = True

        if not allow_fire or not ready:
            safe = max(abs(trigger_value), abs(threshold), 1.0)
            return -safe
        return trigger_value


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
    diagnostics: bool,
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
        elif diagnostics:
            logger.info(
                "delayed_event_fire t=%.6g event=%s trigger=%.6g delay=%.3e",
                time_point,
                scheduled_event.entry.name,
                scheduled_event.trigger_time,
                scheduled_event.delay,
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


def _reconcile_with_time(
    reconcile: Callable[[np.ndarray, float], Dict[str, float]],
    vec: np.ndarray,
    time_point: float,
) -> Dict[str, float]:
    """Helper to call reconcile(vec) or reconcile(vec, t) depending on signature."""
    try:
        return reconcile(vec, time_point)
    except TypeError:
        return reconcile(vec)  # type: ignore[misc]


def _kick_off_t0(
    *,
    model: FrozenModel,
    state: np.ndarray,
    context: Dict[str, float],
    t0: float,
    dt: float,
    reconcile: Callable[[np.ndarray, float], Dict[str, float]],
    method: str = "heun",
    splits: int = 1,
    diagnostics: bool = False,
) -> Tuple[np.ndarray, Dict[str, float], float]:
    """Deterministic explicit advance away from t=0 before running solve_ivp."""

    now = float(t0)
    vec = np.asarray(state, dtype=float, copy=True)
    if dt <= 0.0:
        return vec, context, now
    sub_steps = max(int(splits), 1)
    h = float(dt) / sub_steps
    for _ in range(sub_steps):
        k1 = np.asarray(model.rhs(now, vec.copy()), dtype=float)
        if not np.all(np.isfinite(k1)):
            raise NumericsError("Warm-start kick encountered non-finite RHS (predictor)")
        if method.lower() == "heun":
            y_pred = vec + h * k1
            k2 = np.asarray(model.rhs(now + h, y_pred.copy()), dtype=float)
            if not np.all(np.isfinite(k2)):
                raise NumericsError("Warm-start kick encountered non-finite RHS (corrector)")
            vec = vec + 0.5 * h * (k1 + k2)
        else:
            vec = vec + h * k1
        now = float(now + h)
        context = reconcile(vec, now)
    if diagnostics:
        logger.info("t0_kick: deterministically advanced to t=%.3e", now)
    return vec, context, now



def warm_start_quarantine(
    *,
    model: FrozenModel,
    solver_config: SolverConfig,
    state: np.ndarray,
    context: Dict[str, float],
    reconcile: Callable[[np.ndarray, float], Dict[str, float]],
    samples: Dict[float, np.ndarray],
    time_key: Callable[[float], float],
    pending_events: List[ScheduledEvent],
    t0_options: T0Options,
    jac_sparsity: Optional[np.ndarray],
    tol_time: float,
    diagnostics: bool,
) -> Tuple[np.ndarray, Dict[str, float], float]:
    """Advance snapshot state away from t=0 while suppressing event roots."""

    bump_eps = max(float(t0_options.bump_eps_days or 0.0), np.finfo(float).eps)
    settle = max(float(getattr(t0_options, "settle_no_event_until_days", 0.0) or 0.0), 0.0)
    quarantine = max(float(getattr(t0_options, "quarantine_days", bump_eps) or 0.0), bump_eps)
    target_window = max(bump_eps, settle, quarantine)
    state = np.asarray(state, dtype=float, copy=True)
    context = dict(context)
    current_time = 0.0

    if diagnostics:
        logger.info(
            "t0_warm_start: bump=%.3e settle=%.3e quarantine=%.3e",
            bump_eps,
            settle,
            quarantine,
        )

    def _ensure_rhs_is_finite(t_val: float, vec: np.ndarray) -> None:
        if not t0_options.nan_guard:
            return
        rhs_val = model.rhs(t_val, vec.copy())
        if not np.all(np.isfinite(rhs_val)):
            raise NumericsError(f"RHS non-finite during warm-start (t={t_val:.3e})")

    _ensure_rhs_is_finite(current_time, state)

    triggered_entries: List[EventEntry] = []
    if t0_options.fire_t0_true_events_once:
        for entry in model.events:
            try:
                is_true = _trigger_bool(entry, context)
            except Exception:
                is_true = False
            if is_true:
                triggered_entries.append(entry)

    if triggered_entries and t0_options.fire_t0_true_events_once:
        fire_time = bump_eps
        current_time = fire_time
        if diagnostics:
            logger.info(
                "t0_warm_start: %d events true at t=0 -> firing once at %.3e",
                len(triggered_entries),
                fire_time,
            )
        for entry in triggered_entries:
            delay_val = 0.0
            if entry.delay_compiled is not None:
                delay_val = float(entry.delay_compiled.evaluate(context))
            if entry.delay_type == "time" and delay_val > tol_time:
                scheduled_time = fire_time + delay_val
                heapq.heappush(
                    pending_events,
                    ScheduledEvent(
                        time=scheduled_time,
                        priority=entry.index,
                        entry=entry,
                        trigger_time=current_time,
                        delay=delay_val,
                        assignments=entry.assignments_text,
                    ),
                )
                continue
            for assignment in entry.assignments:
                value = assignment.compiled.evaluate(context) if assignment.compiled is not None else 0.0
                model._apply_target_value(assignment.target, float(value), context, state)
        context = _reconcile_with_time(reconcile, state, current_time)
        samples[time_key(current_time)] = state.copy()

    state, context, current_time = _kick_off_t0(
        model=model,
        state=state,
        context=context,
        t0=current_time,
        dt=bump_eps,
        reconcile=lambda vec, t: _reconcile_with_time(reconcile, vec, t),
        method="heun",
        splits=1,
        diagnostics=diagnostics,
    )
    samples[time_key(current_time)] = state.copy()

    max_trials = max(int(t0_options.tiny_segment_trials or 30), 1)
    max_wall = float(t0_options.max_wall_seconds or 0.0)
    start_wall = time.perf_counter()
    failures = 0

    method_sequence: List[str] = []
    preferred = str(solver_config.method or "").upper()
    if preferred and preferred not in {"RADAU", "BDF"}:
        method_sequence.append(preferred)
    method_sequence.extend(["RADAU", "BDF", "RK23"])
    seen_methods: set = set()
    method_sequence = [
        method
        for method in method_sequence
        if not (method.upper() in seen_methods or seen_methods.add(method.upper()))
    ]

    while current_time < target_window - tol_time:
        if max_wall > 0.0 and (time.perf_counter() - start_wall) > max_wall:
            raise NumericsError(
                f"Warm-start exceeded wall-clock budget (> {max_wall}s) at t={current_time:.3e}"
            )
        if failures >= max_trials:
            raise NumericsError(
                f"Warm-start failed to advance after {max_trials} retries (span≈{target_window-current_time:.3e})"
            )

        span = min(
            target_window - current_time,
            max(t0_options.max_step_cap_days or (target_window - current_time), bump_eps),
        )
        first_step = min(span, t0_options.first_step_cap_days or span)
        max_step = min(span, t0_options.max_step_cap_days or span)

        if diagnostics:
            logger.info(
                "t0_warm_try idx=%d span=%.3e first=%.3e max=%.3e t=%.3e->%.3e",
                failures + 1,
                span,
                first_step,
                max_step,
                current_time,
                current_time + span,
            )

        success = False
        last_reason = ""
        for method_name in method_sequence:
            try:
                sol = solve_ivp(
                    model.rhs,
                    (current_time, current_time + span),
                    state.copy(),
                    method=method_name,
                    rtol=solver_config.rtol,
                    atol=solver_config.atol,
                    first_step=first_step,
                    max_step=max_step,
                    dense_output=False,
                    events=None,
                    jac_sparsity=jac_sparsity if method_name.upper() in {"RADAU", "BDF"} else None,
                    vectorized=False,
                )
            except ValueError as exc:
                last_reason = str(exc)
                if diagnostics:
                    logger.info(
                        "t0_warm_try idx=%d method=%s ValueError: %s",
                        failures + 1,
                        method_name,
                        exc,
                    )
                continue

            if sol.success and sol.y.size and np.isfinite(sol.y).all():
                state = np.asarray(sol.y[:, -1], dtype=float)
                current_time = float(sol.t[-1])
                _ensure_rhs_is_finite(current_time, state)
                context = _reconcile_with_time(reconcile, state, current_time)
                samples[time_key(current_time)] = state.copy()
                failures = 0
                success = True
                if diagnostics:
                    logger.info(
                        "t0_warm_try idx=%d succeeded with method=%s at t=%.3e",
                        failures,
                        method_name,
                        current_time,
                    )
                break
            last_reason = sol.message or "solver failed"
            if diagnostics:
                logger.info(
                    "t0_warm_try idx=%d method=%s failed: %s",
                    failures + 1,
                    method_name,
                    last_reason,
                )

        if success:
            continue

        failures += 1

    return state, context, current_time

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
    reconcile: Callable[[np.ndarray, float], Dict[str, float]],
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
    diagnostics: bool = False,
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
    stagnant_cycles = 0
    event_positions = {spec.entry.index: idx for idx, spec in enumerate(trigger_specs)}
    rearm_below = options.rearm_below if options.rearm_below is not None else -10.0 * abs(options.rearm_eps or 0.0)
    freeze_window = max(options.freeze_window_days, options.settle_no_event_until_days)
    debounce_config = DebounceOptions(
        refractory_dt=options.refractory_dt,
        rearm_requires_false=True,
        rearm_eps=options.rearm_eps,
        rearm_below=rearm_below,
        debounce_window_days=options.debounce_window_days,
        freeze_window_days=freeze_window,
        min_quiet_days=max(options.debounce_window_days, freeze_window),
    )
    debouncer = EventDebouncer(initial_trigger_bools, debounce_config, start_time=start_time)
    limit_first_step = True
    limit_max_step = True

    trigger_functions: List[Callable[[float, np.ndarray], float]] = []
    for spec in trigger_specs:
        pos = event_positions[spec.entry.index]
        original_fn = spec.fn

        def wrapped(t: float, y: np.ndarray, _orig=original_fn, _pos=pos):
            raw_value = float(_orig(t, y))
            return debouncer.observe(_pos, float(t), raw_value, allow_fire=True)

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
            debouncer.observe(pos, time_point, trigger_val, allow_fire=False)

    def reconcile_with_refresh(vec: np.ndarray, t_point: float) -> Dict[str, float]:
        ctx_local = reconcile(vec, t_point)
        refresh_debouncer(t_point, ctx_local)
        return ctx_local

    context = reconcile_with_refresh(state, start_time)
    current_time = start_time

    while segment_index < len(segment_points) - 1 and current_time < stop_time - tol_time:
        cycle_start_time = current_time
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

            if diagnostics:
                logger.info(
                    "segment_solve t=[%.6g, %.6g] pending=%d doses_left=%d",
                    current_time,
                    integration_stop,
                    len(pending_events),
                    len(scheduled_doses) - dose_index,
                )
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
                if diagnostics:
                    logger.warning("segment solver failure at t=%.6g: %s", current_time, sol.message)
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
                if diagnostics:
                    logger.info(
                        "event_fire t=%.6g entries=%s",
                        event_time,
                        ",".join(entry.name for entry in same_time_entries),
                    )
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
                diagnostics=diagnostics,
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
            diagnostics=diagnostics,
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
        if abs(current_time - cycle_start_time) <= tol_time:
            stagnant_cycles += 1
            if stagnant_cycles >= 3:
                if diagnostics:
                    logger.warning("time stagnant near %.6g; forcing nextafter bump", current_time)
                current_time = float(np.nextafter(current_time, np.inf))
                stagnant_cycles = 0
        else:
            stagnant_cycles = 0

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
