import pytest

from src.offline.modules import (
    MODULE_BLOCK_REGISTRY,
    bind_module_blocks,
    disable_repeated_assignments,
    pd1_bridge_block,
    resolve_module_blocks,
    tumour_geometry_block,
)
from src.offline.modules.switches import (
    AVOGADRO,
    DEFAULT_SYN_DEPTH_UM,
    LITERS_PER_CUBIC_MICROMETER,
)
from src.offline.snapshot import load_frozen_model

SURFACE_SCALE = AVOGADRO * LITERS_PER_CUBIC_MICROMETER * DEFAULT_SYN_DEPTH_UM


def test_pd1_bridge_block_combines_compartments():
    model = load_frozen_model("example2")
    block = pd1_bridge_block(model)
    context = {
        "V_T.nivolumab": 2.0,
        "V_P.nivolumab": 1.5,
        "V_C.nivolumab": 0.5,
        "V_LN.nivolumab": 0.25,
        "gamma_T_nivolumab": 0.7,
        "gamma_P_nivolumab": 0.2,
        "gamma_C_nivolumab": 0.05,
        "gamma_LN_nivolumab": 0.1,
    }
    block(context)
    expected_conc = (
        context["gamma_T_nivolumab"] * context["V_T.nivolumab"]
        + context["gamma_P_nivolumab"] * context["V_P.nivolumab"]
        + context["gamma_C_nivolumab"] * context["V_C.nivolumab"]
        + context["gamma_LN_nivolumab"] * context["V_LN.nivolumab"]
    )
    assert context["aPD1_concentration_molar"] == pytest.approx(expected_conc)
    assert context["aPD1"] == pytest.approx(expected_conc * SURFACE_SCALE)


def test_tumour_geometry_block_updates_volume_and_density():
    model = load_frozen_model("example2")
    block = tumour_geometry_block(model)
    context = {
        "V_T": 0.0,
        "cell_to_volume_factor_l": 2e-12,
        "C1": 2.0,
        "C_x": 1.0,
        "V_T.T1": 4.0,
    }
    block(context)
    expected_volume = (context["C1"] + context["C_x"]) * context["cell_to_volume_factor_l"]
    assert context["tumour_volume_l"] == expected_volume
    assert context["tumour_diameter_cm"] > 0.0
    expected_density = context["V_T.T1"] / (expected_volume * 1e6)
    assert context["tcell_density_per_ul"] == expected_density


def test_tumour_geometry_block_derives_volume_from_single_cell_size():
    model = load_frozen_model("example2")
    block = tumour_geometry_block(model)
    context = {
        "V_T": 0.0,
        "vol_cell": 5.0,  # cubic micrometers
        "V_T.C1": 1.5,
        "V_T.C_x": 0.5,
        "V_T.T1": 3.0,
    }
    block(context)
    expected_volume = (context["V_T.C1"] + context["V_T.C_x"]) * 5.0e-15
    assert context["V_T"] == pytest.approx(expected_volume)
    assert context["tumour_volume_l"] == pytest.approx(expected_volume)
    assert context["tcell_density_per_ul"] == pytest.approx(context["V_T.T1"] / (expected_volume * 1e6))


def test_resolve_module_blocks_includes_disable_targets():
    model = load_frozen_model("example2")
    pre_blocks, post_blocks, disable, resolved_names = resolve_module_blocks(
        model, ["pd1_bridge_block"]
    )
    assert resolved_names == ["pd1_bridge_block"]
    assert disable == ["aPD1"]
    assert len(pre_blocks) == 1
    assert len(post_blocks) == 0
    assert "pd1_bridge_block" in MODULE_BLOCK_REGISTRY


def test_resolve_module_blocks_accepts_aliases():
    model = load_frozen_model("example2")
    pre_blocks, post_blocks, disable, resolved_names = resolve_module_blocks(
        model, ["pd1_bridge", "tumour_geometry", "pd1_bridge_block"]
    )
    assert resolved_names == ["pd1_bridge_block", "tumour_geometry_block"]
    assert "aPD1" in disable
    assert len(pre_blocks) == 1
    assert len(post_blocks) == 1


def test_bind_module_blocks_patches_repeated_assignments():
    model = load_frozen_model("example2")
    pre_blocks, post_blocks, disable_targets, _ = resolve_module_blocks(model, ["pd1_bridge_block"])
    disable_repeated_assignments(model, disable_targets)
    bind_module_blocks(model, pre_blocks, post_blocks)

    context = model.build_context_from_state(model.initial_state())
    context.update(
        {
            "V_T.nivolumab": 1.0,
            "V_P.nivolumab": 0.5,
            "gamma_T_nivolumab": 0.8,
            "gamma_P_nivolumab": 0.2,
            "V_C.nivolumab": 0.0,
            "V_LN.nivolumab": 0.0,
            "gamma_C_nivolumab": 0.0,
            "gamma_LN_nivolumab": 0.0,
        }
    )

    model.evaluate_repeated_assignments(context)
    expected_conc = context["gamma_T_nivolumab"] * context["V_T.nivolumab"] + context["gamma_P_nivolumab"] * context["V_P.nivolumab"]
    assert context["aPD1_concentration_molar"] == pytest.approx(expected_conc)
    assert context["aPD1"] == pytest.approx(expected_conc * SURFACE_SCALE)
