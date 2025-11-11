from __future__ import annotations

import numpy as np

from src.offline.simulation import simulate_frozen_model


def test_simulation_records_all_samples_with_warm_start_retry() -> None:
    result = simulate_frozen_model(
        "example1",
        days=0.5,
        therapy="anti_pd1",
        ic_mode="snapshot",
        sample_interval_hours=6.0,
        module_blocks=["pd1_bridge_block", "tumour_geometry_block"],
    )
    expected = np.array([0.0, 0.25, 0.5])
    assert np.allclose(result.time_days, expected)


def test_capture_contexts_returns_debug_payload() -> None:
    result = simulate_frozen_model(
        "example1",
        days=0.5,
        therapy="anti_pd1",
        ic_mode="snapshot",
        sample_interval_hours=6.0,
        module_blocks=["pd1_bridge_block", "tumour_geometry_block"],
        capture_contexts=True,
    )
    assert result.raw_states is not None
    assert result.raw_states.shape[0] == result.time_days.size
    assert result.raw_contexts
    first_ctx = result.raw_contexts[0]
    assert "C1" in first_ctx

