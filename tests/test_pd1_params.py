from __future__ import annotations

from pathlib import Path

import pytest

from src.offline.modules.pd1_params import load_pd1_parameters_from_file
from src.offline.snapshot import load_frozen_model


def test_pd1_parameter_conversion_matches_snapshot_reaction() -> None:
    pd1_params = load_pd1_parameters_from_file(Path("parameters/example1_parameters.json"))
    model = load_frozen_model("artifacts/matlab_frozen_model/example1")
    state = model.initial_state().astype(float)
    model.apply_initial_assignments_to_state(state)
    context = model.build_context_from_state(state.copy())
    model.evaluate_repeated_assignments(context)

    pd1 = float(context["PD1"])
    pdl1 = float(context["PDL1"])
    complex_density = float(context["PD1_PDL1"])
    syn_area = float(context["syn_T1_C1"])

    matlab_kon = float(model.parameters["kon_PD1_PDL1"])
    matlab_koff = float(model.parameters["koff_PD1_PDL1"])

    matlab_flux = (matlab_kon * pd1 * pdl1 - matlab_koff * complex_density) * syn_area
    python_flux = (
        pd1_params.kon_pd1_pdl1 * pd1 * pdl1 - pd1_params.koff_pd1_pdl1 * complex_density
    ) * syn_area

    assert python_flux == pytest.approx(matlab_flux, rel=1e-6)
