from typing import List

import numpy as np

import pandas as pd

from src.offline.frozen_model import EVENT_LOG_FIELDS, _parse_trigger, load_frozen_model, simulate_frozen_model


def test_parse_trigger_direction_downcross():
    expression, direction = _parse_trigger("signal < threshold")
    assert direction == -1.0
    ctx = {"signal": 0.9, "threshold": 1.0}
    value = eval(expression, {}, ctx)  # type: ignore[eval-used]
    assert value < 0.0


def test_event_trigger_direction_with_jitter():
    model = load_frozen_model("example1")
    entry = model.events[0]
    context = {
        "V_T.C1": 0.9 - 1e-9,
        "cell": 1.0,
    }
    # The trigger expression is (V_T.C1)-(0.9*cell); small jitter keeps sign consistent.
    value = entry.trigger_compiled.evaluate(context)
    assert value < 0.0
    flip_context = {
        "V_T.C1": 0.9 + 1e-9,
        "cell": 1.0,
    }
    flipped = entry.trigger_compiled.evaluate(flip_context)
    assert flipped > 0.0


def test_simulation_time_ordered_events_is_stable():
    result = simulate_frozen_model("example2", days=5, therapy="anti_pd1")
    assert np.all(result.t_cells >= 0)


def test_event_log_schema_empty_when_no_events():
    log: List[dict[str, float]] = []
    simulate_frozen_model("example1", days=1, therapy="anti_pd1", event_log=log)
    assert log == []
    frame = pd.DataFrame(log, columns=EVENT_LOG_FIELDS)
    assert list(frame.columns) == list(EVENT_LOG_FIELDS)


def test_event_suite_logs_immediate_and_delayed():
    log: List[dict[str, float]] = []
    simulate_frozen_model("event_suite", days=10.0, therapy="none", event_log=log)
    types = {entry["type"] for entry in log}
    assert types == {"immediate", "delayed"}
    delays = {round(entry["delay"], 9) for entry in log}
    assert 0.0 in delays and 1.0 in delays
