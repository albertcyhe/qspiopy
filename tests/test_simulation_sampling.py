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

