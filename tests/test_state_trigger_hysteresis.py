from __future__ import annotations

import numpy as np

from src.offline.segment_integrator import SolverConfig, T0Options, run_segmented_integration
from src.offline.entities import EventEntry, EventAssignment, CompiledExpression
from src.offline.simulation import _tkey


def _compiled(tokens, func):
    return CompiledExpression(tokens=tuple(tokens), func=func)


class DummyModel:
    def __init__(self, rhs_func):
        self.rhs_func = rhs_func
        self.events = []
        self.fire_count = 0
        self.time_unit = "day"

    def rhs(self, t, y):
        return np.array([self.rhs_func(t, y[0])], dtype=float)

    def build_context_from_state(self, state):
        return {"x": float(state[0])}

    def evaluate_repeated_assignments(self, ctx):
        return

    def apply_algebraic_rules(self, ctx, vec, mutate=False):  # pragma: no cover
        return

    def sync_state_from_context(self, ctx, state):
        state[0] = float(ctx["x"])

    def _apply_target_value(self, target, value, context, state):
        if target == "x":
            state[0] = float(value)
            context["x"] = float(value)
            self.fire_count += 1
        else:
            context[target] = value

    def jacobian_sparsity(self):  # pragma: no cover
        return None


def _make_event_entry(index: int, threshold: float = 0.0) -> EventEntry:
    trigger = _compiled(("x",), lambda x, thr=threshold: x - thr)
    bool_trigger = _compiled(("x",), lambda x, thr=threshold: 1.0 if x - thr >= 0.0 else 0.0)
    assign = EventAssignment(
        target="x",
        expression="x+1",
        compiled=_compiled(("x",), lambda x: x + 1.0),
    )
    return EventEntry(
        index=index,
        name=f"event{index}",
        trigger_expression="x-thr",
        trigger_compiled=trigger,
        trigger_boolean_compiled=bool_trigger,
        direction=1.0,
        delay_expression="0",
        delay_compiled=None,
        delay_type="time",
        assignments=(assign,),
    )


def _run_model(model: DummyModel, event_entry: EventEntry, state0: float, stop_time: float) -> tuple[float, int]:
    trigger_spec = type("Spec", (), {})()
    trigger_spec.entry = event_entry
    trigger_spec.fn = lambda t, y: y[0]
    trigger_spec.fn.direction = 1.0
    trigger_spec.fn.terminal = False
    trigger_spec.state = {"armed": True}

    solver_cfg = SolverConfig(method="BDF", rtol=1e-6, atol=1e-9, max_step=np.inf, seed=None)
    sample_times = np.linspace(0.0, stop_time, 3)
    state = np.array([state0], dtype=float)
    context = model.build_context_from_state(state.copy())
    model.sync_state_from_context(context, state)
    samples = {}

    def record_samples(*args, **kwargs):
        return

    def record_dose(*args, **kwargs):
        return

    def reconcile(vec):
        ctx_local = model.build_context_from_state(vec.copy())
        model.sync_state_from_context(ctx_local, vec)
        return ctx_local

    initial_bool = [bool(event_entry.trigger_boolean_compiled.evaluate(context))]
    t0_opts = T0Options(
        bump_eps_days=0.0,
        debounce_window_days=0.0,
        nan_guard=False,
    )

    final_state, _, _, _ = run_segmented_integration(
        model=model,
        solver_config=solver_cfg,
        scheduled_doses=[],
        pending_events=[],
        trigger_specs=[trigger_spec],
        initial_trigger_bools=initial_bool,
        sample_times=sample_times,
        samples=samples,
        state=state,
        context=context,
        start_time=0.0,
        stop_time=stop_time,
        dose_index=0,
        tol_time=1e-9,
        record_solution_samples=record_samples,
        record_dose_audit=record_dose,
        reconcile=reconcile,
        time_key=_tkey,
        event_log=None,
        jac_sparsity=model.jacobian_sparsity(),
        record_state=None,
        t0_options=t0_opts,
    )
    return float(final_state[0]), model.fire_count


def test_event_prearmed_does_not_fire():
    model = DummyModel(rhs_func=lambda t, x: 0.0)
    event = _make_event_entry(index=0, threshold=0.0)
    final_x, fire_count = _run_model(model, event, state0=0.0, stop_time=0.1)
    assert fire_count == 0
    assert final_x == 0.0


def test_event_rearms_after_crossing():
    model = DummyModel(rhs_func=lambda t, x: 1.0)
    event = _make_event_entry(index=0, threshold=0.0)
    final_x, fire_count = _run_model(model, event, state0=-0.5, stop_time=1.5)
    assert fire_count == 1
    assert final_x > 0.0
