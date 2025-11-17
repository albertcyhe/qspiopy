"""PD-1 parameter loading helpers with explicit unit conversions."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, MutableMapping, Optional, Tuple

from ..units import convert_area, convert_length

SECONDS_PER_DAY = 86400.0


@dataclass(frozen=True)
class PD1Params:
    kon_pd1_pdl1: float
    kon_pd1_pdl2: float
    kon_pd1_ab: float
    koff_pd1_pdl1: float
    koff_pd1_pdl2: float
    koff_pd1_ab: float
    chi_pd1: float
    internalisation_per_day: float
    smoothing_tau_days: float
    pd1_50_density: float
    hill_coefficient: float
    gamma_c_nivolumab: float
    synapse_area_um2: float
    tcell_area_um2: float
    cell_area_um2: float
    total_pd1_molecules: float
    total_pdl1_molecules: float
    total_pdl2_molecules: float
    max_step_days: float
    solver_rtol: float
    solver_atol: float

    def pd1_surface_density(self) -> float:
        area = max(self.tcell_area_um2, 1e-6)
        return self.total_pd1_molecules / area

    def pdl1_surface_density(self) -> float:
        area = max(self.cell_area_um2, 1e-6)
        return self.total_pdl1_molecules / area

    def pdl2_surface_density(self) -> float:
        area = max(self.cell_area_um2, 1e-6)
        return self.total_pdl2_molecules / area


def _read_parameter_entries(path: Path) -> Dict[str, Dict[str, object]]:
    data = json.loads(path.read_text())
    entries: Dict[str, Dict[str, object]] = {}
    for entry in data:
        name = entry.get("name")
        if name:
            entries[str(name)] = entry
    return entries


def _entry_value(entries: Mapping[str, Dict[str, object]], name: str) -> float:
    entry = entries.get(name)
    if entry is None or entry.get("value") is None:
        raise KeyError(f"Parameter '{name}' missing from parameter file")
    return float(entry["value"])


def _entry_units(entries: Mapping[str, Dict[str, object]], name: str) -> str:
    entry = entries.get(name)
    if entry is None:
        return ""
    units = entry.get("units")
    return str(units or "").lower()


def _length_to_micrometer(value: float, unit: str) -> float:
    meters = convert_length(value, unit or "meter")
    return meters / 1e-6


def _area_to_micrometer2(value: float, unit: str) -> float:
    area_m2 = convert_area(value, unit or "meter^2")
    return area_m2 / 1e-12


def _optional_value(entries: Mapping[str, Dict[str, object]], name: str, default: float = 0.0) -> float:
    entry = entries.get(name)
    if entry is None or entry.get("value") is None:
        return default
    return float(entry["value"])


def _surface_area(entries: Mapping[str, Dict[str, object]], area_name: str, diameter_name: str) -> float:
    entry = entries.get(area_name)
    if entry and entry.get("value") is not None:
        return _area_to_micrometer2(float(entry["value"]), _entry_units(entries, area_name))
    diameter = _entry_value(entries, diameter_name)
    units = _entry_units(entries, diameter_name)
    diameter_um = _length_to_micrometer(diameter, units)
    return 4.0 * math.pi * (diameter_um / 2.0) ** 2


def _synapse_depth_um(entries: Mapping[str, Dict[str, object]]) -> float:
    depth = _entry_value(entries, "d_syn")
    units = _entry_units(entries, "d_syn")
    return _length_to_micrometer(depth, units)


def _kd_to_molar(value: float, unit: str) -> float:
    norm = (unit or "").lower()
    if "nanomolar" in norm:
        return value * 1e-9
    if "micromolar" in norm:
        return value * 1e-6
    return value


def _kd_to_micromolar(value: float, unit: str) -> float:
    norm = (unit or "").lower()
    if "nanomolar" in norm:
        return value * 1e-3
    if "micromolar" in norm:
        return value
    return value


def load_pd1_parameters_from_file(path: Path) -> PD1Params:
    return load_pd1_parameters_from_entries(_read_parameter_entries(path))


def load_pd1_parameters_from_entries(entries: Mapping[str, Dict[str, object]]) -> PD1Params:
    d_syn_um = _synapse_depth_um(entries)
    d_syn_nm = d_syn_um * 1000.0
    synapse_area_um2 = _area_to_micrometer2(_entry_value(entries, "A_syn"), _entry_units(entries, "A_syn"))
    gamma_c = _optional_value(entries, "gamma_C_nivolumab", 1.0)
    internalisation = _optional_value(entries, "pd1_occ_internalization_per_day", 0.0)
    smoothing_tau = _optional_value(entries, "pd1_whitebox_tau_days", 0.0)
    pd1_50 = _optional_value(entries, "PD1_50", 1.0)
    hill_coeff = _optional_value(entries, "n_PD1", 1.0)
    max_step = _optional_value(entries, "pd1_whitebox_max_step_days", 1e-4)
    solver_rtol = _optional_value(entries, "pd1_whitebox_rtol", 1e-6)
    solver_atol = _optional_value(entries, "pd1_whitebox_atol", 1e-9)

    k_pd1_pdl1 = _entry_value(entries, "k_PD1_PDL1")
    k_pd1_pdl2 = _entry_value(entries, "k_PD1_PDL2")
    kd_pd1_pdl1 = _entry_value(entries, "kd_PD1_PDL1")
    kd_pd1_pdl2 = _entry_value(entries, "kd_PD1_PDL2")
    kd_pd1_pdl1_u = _entry_units(entries, "kd_PD1_PDL1")
    kd_pd1_pdl2_u = _entry_units(entries, "kd_PD1_PDL2")

    kon_pd1_pdl1 = k_pd1_pdl1 * SECONDS_PER_DAY / (max(d_syn_nm, 1e-9) * 1e3)
    kon_pd1_pdl2 = k_pd1_pdl2 * SECONDS_PER_DAY / (max(d_syn_nm, 1e-9) * 1e3)

    koff_pd1_pdl1 = k_pd1_pdl1 * kd_pd1_pdl1 * SECONDS_PER_DAY
    koff_pd1_pdl2 = k_pd1_pdl2 * kd_pd1_pdl2 * SECONDS_PER_DAY

    k_pd1_ab = _entry_value(entries, "kon_PD1_aPD1")
    kd_pd1_ab = _entry_value(entries, "kd_PD1_aPD1")
    kd_pd1_ab_u = _entry_units(entries, "kd_PD1_aPD1")
    kon_pd1_ab = k_pd1_ab * SECONDS_PER_DAY
    koff_pd1_ab = kon_pd1_ab * _kd_to_molar(kd_pd1_ab, kd_pd1_ab_u)

    chi_pd1 = _entry_value(entries, "chi_PD1") / max(d_syn_nm, 1e-9)

    total_pd1 = _entry_value(entries, "T_PD1_total")
    total_pdl1 = _entry_value(entries, "C_PDL1_total")
    total_pdl2 = _entry_value(entries, "C_PDL2_total")

    tcell_area = _surface_area(entries, "A_Tcell", "D_Tcell")
    cell_area = _surface_area(entries, "A_cell", "D_cell")

    return PD1Params(
        kon_pd1_pdl1=kon_pd1_pdl1,
        kon_pd1_pdl2=kon_pd1_pdl2,
        kon_pd1_ab=kon_pd1_ab,
        koff_pd1_pdl1=koff_pd1_pdl1,
        koff_pd1_pdl2=koff_pd1_pdl2,
        koff_pd1_ab=koff_pd1_ab,
        chi_pd1=chi_pd1,
        internalisation_per_day=float(internalisation),
        smoothing_tau_days=float(smoothing_tau),
        pd1_50_density=float(pd1_50),
        hill_coefficient=float(hill_coeff),
        gamma_c_nivolumab=float(gamma_c),
        synapse_area_um2=float(synapse_area_um2),
        tcell_area_um2=float(tcell_area),
        cell_area_um2=float(cell_area),
        total_pd1_molecules=float(total_pd1),
        total_pdl1_molecules=float(total_pdl1),
        total_pdl2_molecules=float(total_pdl2),
        max_step_days=max(float(max_step), 1e-8),
        solver_rtol=max(float(solver_rtol), 1e-12),
        solver_atol=max(float(solver_atol), 1e-18),
    )


def pd1_params_from_snapshot(parameters: Mapping[str, float]) -> PD1Params:
    return PD1Params(
        kon_pd1_pdl1=float(parameters.get("kon_PD1_PDL1", 0.0)),
        kon_pd1_pdl2=float(parameters.get("kon_PD1_PDL2", 0.0)),
        kon_pd1_ab=float(parameters.get("kon_PD1_aPD1", 0.0)),
        koff_pd1_pdl1=float(parameters.get("koff_PD1_PDL1", 0.0)),
        koff_pd1_pdl2=float(parameters.get("koff_PD1_PDL2", 0.0)),
        koff_pd1_ab=float(parameters.get("koff_PD1_aPD1", 0.0)),
        chi_pd1=float(parameters.get("Chi_PD1", 0.0)),
        internalisation_per_day=float(parameters.get("pd1_occ_internalization_per_day", 0.0)),
        smoothing_tau_days=float(parameters.get("pd1_whitebox_tau_days", 0.0)),
        pd1_50_density=float(parameters.get("pd1_whitebox_pd1_50_density", parameters.get("PD1_50", 1.0))),
        hill_coefficient=float(parameters.get("n_PD1", 1.0)),
        gamma_c_nivolumab=float(parameters.get("gamma_C_nivolumab", 1.0)),
        synapse_area_um2=float(parameters.get("A_syn", 1.0)),
        tcell_area_um2=float(parameters.get("A_Tcell", 1.0)),
        cell_area_um2=float(parameters.get("A_cell", 1.0)),
        total_pd1_molecules=float(parameters.get("T_PD1_total", 0.0)),
        total_pdl1_molecules=float(parameters.get("C_PDL1_total", 0.0)),
        total_pdl2_molecules=float(parameters.get("C_PDL2_total", 0.0)),
        max_step_days=float(parameters.get("pd1_whitebox_max_step_days", 1e-4) or 1e-4),
        solver_rtol=float(parameters.get("pd1_whitebox_rtol", 1e-6) or 1e-6),
        solver_atol=float(parameters.get("pd1_whitebox_atol", 1e-9) or 1e-9),
    )


__all__ = [
    "PD1Params",
    "load_pd1_parameters_from_entries",
    "load_pd1_parameters_from_file",
    "pd1_params_from_snapshot",
]
