from __future__ import annotations

import math

import pytest

from src.offline.units import (
    convert_parameter_value,
    flow_to_l_per_day,
    normalise_dose_to_species,
)


def test_flow_to_l_per_day_converts_from_ml_hour() -> None:
    assert flow_to_l_per_day(1.0, "ml/hour") == pytest.approx(24e-3)
    assert flow_to_l_per_day(2.5, "l/day") == pytest.approx(2.5)


def test_convert_parameter_value_kon_legacy_factor() -> None:
    converted = convert_parameter_value(2.0, "1/(micromolarity*nanometer*second)")
    assert converted == pytest.approx(2.0 * 9.8412890625)


def test_normalise_dose_to_species_handles_concentration_targets() -> None:
    delta, amount_mol = normalise_dose_to_species(
        amount=200.0,
        amount_unit="milligram",
        species_units="mole/liter",
        species_dimension="concentration",
        compartment_volume_l=5.0,
        molecular_weight_g_per_mol=1.436e5,
    )
    expected_moles = 0.2 / 1.436e5
    assert amount_mol == pytest.approx(expected_moles)
    assert delta == pytest.approx(expected_moles / 5.0)


def test_normalise_dose_to_species_prefers_units_over_dimension() -> None:
    delta, amount_mol = normalise_dose_to_species(
        amount=200.0,
        amount_unit="milligram",
        species_units="molarity",
        species_dimension="amount",
        compartment_volume_l=5.0,
        molecular_weight_g_per_mol=1.436e5,
    )
    expected_moles = 0.2 / 1.436e5
    assert delta == pytest.approx(expected_moles / 5.0)
    assert amount_mol == pytest.approx(expected_moles)
