import pytest

from src.offline.modules import pd1_occupancy_filter_block
from src.offline.snapshot import load_frozen_model


def _build_pd1_block(**overrides):
    model = load_frozen_model("example2")
    model.parameters = dict(model.parameters)
    model.parameters.update(overrides)
    return pd1_occupancy_filter_block(model)


def test_pd1_occupancy_filter_filters_signal():
    block = _build_pd1_block()

    ctx0 = {
        "time_days": 0.0,
        "aPD1_surface_molecules_per_um2": 0.0,
        "tumour_volume_l": 0.01,
    }
    block(ctx0)
    assert ctx0["pd1_occupancy"] == pytest.approx(0.0)
    assert ctx0["H_PD1_C1"] == pytest.approx(0.0)

    ctx1 = {
        "time_days": 1.0,
        "aPD1_surface_molecules_per_um2": 100.0,
        "tumour_volume_l": 0.01,
    }
    block(ctx1)
    assert 0.0 < ctx1["pd1_occupancy"] < 1.0

    ctx2 = {
        "time_days": 10.0,
        "aPD1_surface_molecules_per_um2": 100.0,
        "tumour_volume_l": 0.01,
    }
    block(ctx2)
    assert ctx2["pd1_occupancy"] >= ctx1["pd1_occupancy"]
    assert ctx2["H_PD1_C1"] == ctx2["pd1_occupancy"]


def test_pd1_occupancy_filter_applies_delay_steps():
    block = _build_pd1_block(pd1_occ_delay_steps=2, pd1_occ_tau1_days=1.0, pd1_occ_tau2_days=1.0)
    outputs = []
    ctx = {"time_days": 0.0, "aPD1_surface_molecules_per_um2": 0.0, "tumour_volume_l": 0.01}
    block(ctx)
    outputs.append(ctx["pd1_occupancy"])
    ctx = {"time_days": 1.0, "aPD1_surface_molecules_per_um2": 200.0, "tumour_volume_l": 0.01}
    block(ctx)
    outputs.append(ctx["pd1_occupancy"])
    ctx = {"time_days": 2.0, "aPD1_surface_molecules_per_um2": 200.0, "tumour_volume_l": 0.01}
    block(ctx)
    outputs.append(ctx["pd1_occupancy"])
    ctx = {"time_days": 3.0, "aPD1_surface_molecules_per_um2": 200.0, "tumour_volume_l": 0.01}
    block(ctx)
    outputs.append(ctx["pd1_occupancy"])
    assert outputs[1] == pytest.approx(0.0)  # still delayed
    assert outputs[-1] > outputs[1]
