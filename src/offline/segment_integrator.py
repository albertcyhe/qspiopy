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

_MICRO_SNAP = 5e-9


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


def _process_pending_events(
    *,
    time_point: float,
    pending_events: List[ScheduledEvent],
    model: FrozenModel,
    context: Dict[str, float],
    state: np.ndarray,
    reconcile: Callable[[np.ndarray], Dict[str, float]],
    samples: Dict[float, np.ndarray],
    time_key: Callable[[float], float],
    event_log: Optional[List[Dict[str, object]]],
    tol_time: float,
) -> None:
    while pending_events and abs(pending_events[0].time - time_point) <= tol_time:
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
        context = reconcile(state)
        samples[time_key(time_point)] = state.copy()


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


def run_segmented_integration(
    *,
    model: FrozenModel,
    solver_config: SolverConfig,
    scheduled_doses: Sequence[ScheduledDose],
    pending_events: List[ScheduledEvent],
    trigger_specs: Sequence[TriggerEventSpec],
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
    current_time = start_time
    segment_index = 0
    trigger_functions = [spec.fn for spec in trigger_specs]
    trigger_entries = [spec.entry for spec in trigger_specs]

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

            sol = solve_ivp(
                model.rhs,
                (current_time, integration_stop),
                state,
                method=solver_config.method,
                rtol=solver_config.rtol,
                atol=solver_config.atol,
                max_step=solver_config.max_step,
                first_step=_initial_step(current_time, integration_stop),
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
                context = reconcile(state)
                same_time_entries = [
                    entry for time_val, entry in triggered_pairs if abs(time_val - event_time) <= tol_time
                ]
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
                        model._apply_target_value(assignment.target, value, context, state)
                    context = reconcile(state)
                samples[time_key(event_time)] = state.copy()
                if record_state:
                    last_entry = same_time_entries[-1] if same_time_entries else triggered_pairs[0][1]
                    record_state(event_time, state.copy(), "post", "event", last_entry.index, last_entry.name, None, None)
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
            context = reconcile(state)
            current_time = integration_stop

            _process_pending_events(
                time_point=current_time,
                pending_events=pending_events,
                model=model,
                context=context,
                state=state,
                reconcile=reconcile,
                samples=samples,
                time_key=time_key,
                event_log=event_log,
                tol_time=tol_time,
            )

            if abs(current_time - target_time) <= tol_time:
                reached_target = True

        current_time = target_time

        _process_pending_events(
            time_point=current_time,
            pending_events=pending_events,
            model=model,
            context=context,
            state=state,
            reconcile=reconcile,
            samples=samples,
            time_key=time_key,
            event_log=event_log,
            tol_time=tol_time,
        )

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
                context = reconcile(state)
                samples[time_key(current_time)] = state.copy()
                record_dose_audit(current_time, scheduled_dose, audit_info)
                if record_state:
                    record_state(current_time, state.copy(), "post", "dose", None, scheduled_dose.dose.name, scheduled_dose.dose.target, scheduled_dose.dose.name)

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
