from __future__ import annotations

from pathlib import Path

import pytest

from src.offline.modules.pd1_params import load_pd1_parameters_from_file


def test_pd1_parameter_conversion_matches_expected_units() -> None:
    pd1_params = load_pd1_parameters_from_file(Path("parameters/example1_parameters.json"))
    # Derived from k_PD1_PDL1 / d_syn (0.175 / 3 nm) and k_PD1_PDL1 * kd_PD1_PDL1.
    assert pd1_params.kon_pd1_pdl1 == pytest.approx(0.05833333333333333, rel=1e-12)
    assert pd1_params.koff_pd1_pdl1 == pytest.approx(1.435, rel=1e-12)
