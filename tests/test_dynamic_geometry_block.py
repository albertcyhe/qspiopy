import math

import pytest

from src.offline.modules import tumour_geometry_dynamic_block
from src.offline.snapshot import load_frozen_model


def _build_block():
    model = load_frozen_model("example2")
    return tumour_geometry_dynamic_block(model)


def _build_block_with_params(**overrides):
    model = load_frozen_model("example2")
    model.parameters = dict(model.parameters)
    model.parameters.update(overrides)
    return tumour_geometry_dynamic_block(model)


def test_dynamic_geometry_block_aggregates_clones_and_updates_density():
    block = _build_block()
    context = {
        "time_days": 0.0,
        "cell_to_volume_factor_l": 2e-12,
        "V_T.C1": 2.0,
        "V_T.C2": 3.0,
        "V_T.T1": 5.0,
        "V_T.C_x1": 1.0,
    }
    block(context)
    v_cell = 2e-12
    live_volume = (2.0 + 3.0) * v_cell
    dead_volume = 1.0 * 1.2 * v_cell
    t_volume = 5.0 * v_cell
    expected_total = live_volume + dead_volume + t_volume
    assert context["tumour_volume_l"] == pytest.approx(expected_total)
    assert context["tumour_diameter_cm"] > 0.0
    expected_density = 5.0 / (expected_total * 1e6)
    assert context["tcell_density_per_ul"] == pytest.approx(expected_density)


def test_dynamic_geometry_block_applies_dead_volume_filter():
    block = _build_block()
    context = {
        "time_days": 0.0,
        "cell_to_volume_factor_l": 1e-12,
        "V_T.C1": 1.0,
        "V_T.C_x1": 2.0,
    }
    block(context)
    first_filtered = context["geom_dead_filtered_volume_l"]
    context.update({"time_days": 1.0, "V_T.C_x1": 0.0})
    block(context)
    second_filtered = context["geom_dead_filtered_volume_l"]
    assert second_filtered < first_filtered
    assert second_filtered > 0.0


def test_dynamic_geometry_block_uses_litre_scaled_vol_cell_when_small():
    block = _build_block()
    context = {
        "time_days": 0.0,
        "vol_cell": 2e-12,  # already in litres
        "V_T.C1": 1.0,
    }
    block(context)
    assert context["geom_live_volume_l"] == pytest.approx(2e-12)


def test_dynamic_geometry_block_relaxes_toward_target_volume():
    block = _build_block_with_params(geom_tau_days=10.0)
    ctx0 = {
        "time_days": 0.0,
        "tumour_volume_l": 1.0,
        "vol_cell": 1e-12,
    }
    block(ctx0)
    initial_volume = ctx0["tumour_volume_l"]
    assert initial_volume == pytest.approx(1.0)

    ctx1 = {
        "time_days": 1.0,
        "tumour_volume_l": 5.0,
        "vol_cell": 1e-12,
    }
    block(ctx1)
    updated_volume = ctx1["tumour_volume_l"]
    assert initial_volume < updated_volume < 5.0
