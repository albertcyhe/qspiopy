"""Unit conversion helpers for frozen SimBiology snapshots."""

from __future__ import annotations

import math
from typing import Literal

TimeUnit = Literal["day"]
VolumeUnit = Literal["liter"]
AmountUnit = Literal["mole"]
ConcentrationUnit = Literal["molarity"]


def _normalize(unit: str) -> str:
    return (unit or "").strip().lower()


_VOLUME_FACTORS = {
    "litre": 1.0,
    "liter": 1.0,
    "milliliter": 1e-3,
    "millilitre": 1e-3,
    "microliter": 1e-6,
    "microlitre": 1e-6,
    "millimeter^3": 1e-6,
    "micrometer^3": 1e-15,
    "micrometer^3/cell": 1e-15,
    "centimeter^3": 1e-3,
}


def convert_volume(value: float, unit: str) -> float:
    norm = _normalize(unit)
    if norm in _VOLUME_FACTORS:
        return value * _VOLUME_FACTORS[norm]
    return value


_CONCENTRATION_FACTORS = {
    "molarity": 1.0,
    "micromolarity": 1e-6,
    "nanomolarity": 1e-9,
    "mole/liter": 1.0,
    "mol/l": 1.0,
}


def convert_concentration(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _CONCENTRATION_FACTORS.get(norm)
    if factor is not None:
        return value * factor
    return value



def convert_kon(value: float, unit: str) -> float:
    norm = _normalize(unit)
    if not norm:
        return value
    if norm == "1/(molarity*second)" or norm == "1/(mole/liter*second)":
        return value * 86400.0
    if norm == "1/(micromolarity*nanometer*second)":
        # uM^-1 nm^-1 s^-1 to L/(mol*day)
        return value * 9.8412890625
    if norm == "1/(micromolarity*second)":
        return value * 86400.0 * 1e6
    return value


_RATE_FACTORS = {
    "1/day": 1.0,
    "1/minute": 1440.0,
    "1/min": 1440.0,
    "1/second": 86400.0,
    "1/sec": 86400.0,
    "1/hour": 24.0,
    "1/h": 24.0,
    "nanomole/cell/hour": 1e-9 * 24.0,
    "1/(centimeter^3*minute)": 1e3 * 1440.0,
    "1/day/milliliter": 1000.0,
}


def convert_rate(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _RATE_FACTORS.get(norm)
    if factor is not None:
        return value * factor
    if norm.endswith("per day"):
        return value
    return value


def _convert_denominator(unit: str) -> float:
    norm = _normalize(unit)
    if not norm or norm in {"cell", "cells", "molecule", "molecules"}:
        return 1.0
    if norm in _VOLUME_FACTORS:
        return convert_volume(1.0, norm)
    if norm in _AREA_FACTORS:
        return convert_area(1.0, norm)
    if norm in _LENGTH_FACTORS:
        return convert_length(1.0, norm)
    if norm == "molarity":
        return convert_concentration(1.0, "molarity")
    return 1.0


def convert_parameter_value(value: float, unit: str) -> float:
    norm = _normalize(unit)
    if not norm or norm == "dimensionless":
        return value
    if norm in _VOLUME_FACTORS:
        return convert_volume(value, unit)
    if norm in _CONCENTRATION_FACTORS:
        return convert_concentration(value, unit)
    if norm.startswith("1/(mol") or norm.startswith("1/(micromolar"):
        return convert_kon(value, unit)
    ratio_prefixes = {
        "cell/day": lambda v: v * 1.0 / convert_parameter_value(1.0, "cell"),
        "molecule/day": lambda v: v * 1.0 / convert_parameter_value(1.0, "molecule"),
    }
    if norm in ratio_prefixes:
        return ratio_prefixes[norm](value)

    if any(norm.startswith(prefix) for prefix in ("1/day/", "1/hour/", "1/minute/", "1/second/")):
        for prefix in ("1/day/", "1/hour/", "1/minute/", "1/second/"):
            if norm.startswith(prefix):
                denom = norm[len(prefix):]
                return convert_rate(value, prefix.rstrip("/")) / _convert_denominator(denom)
    if norm.startswith("1/") and norm[2:] in _LENGTH_FACTORS:
        return value / convert_length(1.0, norm[2:])
    if norm.startswith("1/(second"):
        return convert_rate(value, "1/second")
    if norm in _RATE_FACTORS or norm.startswith("1/"):
        return convert_rate(value, unit)
    if norm == "cell/milliliter":
        return value * 1000.0
    return value


_MOLE_FACTORS = {
    "mole": 1.0,
    "mol": 1.0,
    "millimole": 1e-3,
    "mmol": 1e-3,
    "micromole": 1e-6,
    "µmol": 1e-6,
}

_MASS_TO_GRAM = {
    "gram": 1.0,
    "g": 1.0,
    "milligram": 1e-3,
    "mg": 1e-3,
    "microgram": 1e-6,
    "ug": 1e-6,
    "µg": 1e-6,
}


def convert_amount(value: float, unit: str, molecular_weight_g_per_mol: float | None = None) -> float:
    """Convert `value` with given unit into moles."""
    norm = _normalize(unit)
    if not norm:
        return value
    if norm in _MOLE_FACTORS:
        return value * _MOLE_FACTORS[norm]
    if norm in _MASS_TO_GRAM:
        if molecular_weight_g_per_mol is None or molecular_weight_g_per_mol == 0.0:
            raise ValueError("Molecular weight (g/mol) required to convert mass to moles")
        grams = value * _MASS_TO_GRAM[norm]
        return grams / molecular_weight_g_per_mol
    return value
_LENGTH_FACTORS = {
    "meter": 1.0,
    "metre": 1.0,
    "centimeter": 1e-2,
    "centimetre": 1e-2,
    "millimeter": 1e-3,
    "millimetre": 1e-3,
    "micrometer": 1e-6,
    "micrometre": 1e-6,
    "nanometer": 1e-9,
    "nanometre": 1e-9,
}


def convert_length(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _LENGTH_FACTORS.get(norm)
    if factor is not None:
        return value * factor
    return value


_AREA_FACTORS = {
    "meter^2": 1.0,
    "metre^2": 1.0,
    "centimeter^2": 1e-4,
    "centimetre^2": 1e-4,
    "millimeter^2": 1e-6,
    "micrometer^2": 1e-12,
    "micrometre^2": 1e-12,
    "nanometer^2": 1e-18,
    "nanometre^2": 1e-18,
}


def convert_area(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _AREA_FACTORS.get(norm)
    if factor is not None:
        return value * factor
    return value
