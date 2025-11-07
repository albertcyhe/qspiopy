"""Unit conversion helpers for frozen SimBiology snapshots."""

from __future__ import annotations

from typing import Optional, Tuple

DAY_PER_SECOND = 86400.0
DAY_PER_MINUTE = 1440.0
DAY_PER_HOUR = 24.0


def _normalize(unit: str) -> str:
    text = (unit or "").strip().lower()
    return text.replace("µ", "u").replace(" ", "")


_VOLUME_FACTORS = {
    "l": 1.0,
    "liter": 1.0,
    "litre": 1.0,
    "milliliter": 1e-3,
    "millilitre": 1e-3,
    "ml": 1e-3,
    "microliter": 1e-6,
    "microlitre": 1e-6,
    "ul": 1e-6,
    "nanoliter": 1e-9,
    "nanolitre": 1e-9,
    "nl": 1e-9,
    "millimeter^3": 1e-6,
    "millimetre^3": 1e-6,
    "mm^3": 1e-6,
    "centimeter^3": 1e-3,
    "centimetre^3": 1e-3,
    "cm^3": 1e-3,
    "micrometer^3": 1e-15,
    "micrometre^3": 1e-15,
    "um^3": 1e-15,
    "micrometer^3/cell": 1e-15,
    "micrometre^3/cell": 1e-15,
    "um^3/cell": 1e-15,
}


def convert_volume(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _VOLUME_FACTORS.get(norm)
    if factor is not None:
        return value * factor
    return float(value)


_CONCENTRATION_FACTORS = {
    "m": 1.0,
    "molar": 1.0,
    "molarity": 1.0,
    "mol/l": 1.0,
    "mole/liter": 1.0,
    "millimolar": 1e-3,
    "mm": 1e-3,
    "micromolar": 1e-6,
    "micromolarity": 1e-6,
    "um": 1e-6,
    "nanomolar": 1e-9,
    "nanomolarity": 1e-9,
    "nm": 1e-9,
}


def convert_concentration(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _CONCENTRATION_FACTORS.get(norm)
    if factor is not None:
        return value * factor
    return float(value)


def time_to_day(value: float, unit: str) -> float:
    norm = _normalize(unit)
    if norm in {"day", "days", "d"}:
        return float(value)
    if norm in {"hour", "hours", "h"}:
        return float(value) / DAY_PER_HOUR
    if norm in {"minute", "minutes", "min"}:
        return float(value) / DAY_PER_MINUTE
    if norm in {"second", "seconds", "sec", "s"}:
        return float(value) / DAY_PER_SECOND
    if norm in {"week", "weeks", "wk"}:
        return float(value) * 7.0
    return float(value)


def rate_to_per_day(value: float, unit: str) -> float:
    norm = _normalize(unit)
    if not norm or norm == "1/day" or norm == "perday" or norm == "day^-1" or norm == "d^-1":
        return float(value)
    if norm in {"1/hour", "perhour", "hour^-1", "h^-1"}:
        return float(value) * DAY_PER_HOUR
    if norm in {"1/minute", "perminute", "minute^-1", "min^-1"}:
        return float(value) * DAY_PER_MINUTE
    if norm in {"1/second", "persecond", "second^-1", "sec^-1", "1/sec", "1/s"}:
        return float(value) * DAY_PER_SECOND
    if norm in {"1/week", "perweek", "week^-1"}:
        return float(value) / 7.0
    return float(value)


_FLOW_FACTORS = {
    "l/day": 1.0,
    "liter/day": 1.0,
    "litre/day": 1.0,
    "l/d": 1.0,
    "ml/day": 1e-3,
    "milliliter/day": 1e-3,
    "millilitre/day": 1e-3,
    "ml/d": 1e-3,
    "l/hour": 24.0,
    "liter/hour": 24.0,
    "litre/hour": 24.0,
    "l/h": 24.0,
    "ml/hour": 24.0 * 1e-3,
    "milliliter/hour": 24.0 * 1e-3,
    "millilitre/hour": 24.0 * 1e-3,
    "ml/h": 24.0 * 1e-3,
    "l/minute": 1440.0,
    "liter/minute": 1440.0,
    "litre/minute": 1440.0,
    "ml/minute": 1440.0 * 1e-3,
}


def flow_to_l_per_day(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _FLOW_FACTORS.get(norm)
    if factor is not None:
        return float(value) * factor
    return float(value)


_SPECIAL_RATE_FACTORS = {
    "nanomole/cell/hour": 1e-9 * DAY_PER_HOUR,
    "nanomole/cell/day": 1e-9,
    "nanomole/cell/minute": 1e-9 * DAY_PER_MINUTE,
    "molecule/micrometer^2/second": DAY_PER_SECOND / (1e-12),
    "molecule/micrometre^2/second": DAY_PER_SECOND / (1e-12),
    "1/(centimeter^3*minute)": 1e3 * DAY_PER_MINUTE,
    "1/(centimetre^3*minute)": 1e3 * DAY_PER_MINUTE,
    "1/(cm^3*minute)": 1e3 * DAY_PER_MINUTE,
}


def convert_rate(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _SPECIAL_RATE_FACTORS.get(norm)
    if factor is not None:
        return float(value) * factor
    if norm.endswith("perday"):
        return float(value)
    if norm.startswith("1/") or norm.endswith("^-1"):
        return rate_to_per_day(value, unit)
    return float(value)


def convert_kon(value: float, unit: str) -> float:
    norm = _normalize(unit)
    if not norm:
        return float(value)
    if norm in {"1/(molarity*second)", "1/(mol/liter*second)", "1/(m*s)", "1/(molar*s)"}:
        return float(value) * DAY_PER_SECOND
    if norm in {
        "1/(micromolarity*nanometer*second)",
        "1/(um*nm*second)",
        "1/(micromolarity*nm*second)",
    }:
        return float(value) * 9.8412890625
    if norm in {"1/(micromolarity*second)", "1/(um*second)"}:
        return float(value) * DAY_PER_SECOND * 1e6
    return float(value)


def _convert_denominator(unit: str) -> float:
    norm = _normalize(unit)
    if "/" in norm:
        result = 1.0
        for part in norm.split("/"):
            if part:
                result *= _convert_denominator(part)
        return result
    if not norm or norm in {"cell", "cells", "molecule", "molecules", "count"}:
        return 1.0
    if norm in {"day", "days", "d"}:
        return 1.0
    if norm in {"hour", "hours", "h", "minute", "minutes", "min", "second", "seconds", "sec", "s", "week", "weeks", "wk"}:
        return time_to_day(1.0, norm)
    if norm in _VOLUME_FACTORS:
        return convert_volume(1.0, norm)
    if norm in _AREA_FACTORS:
        return convert_area(1.0, norm)
    if norm in _LENGTH_FACTORS:
        return convert_length(1.0, norm)
    if norm in _CONCENTRATION_FACTORS:
        return convert_concentration(1.0, norm)
    if norm in _FLOW_FACTORS:
        return flow_to_l_per_day(1.0, norm)
    if norm in _SPECIAL_RATE_FACTORS:
        return _SPECIAL_RATE_FACTORS[norm]
    return 1.0


def convert_parameter_value(value: float, unit: str) -> float:
    norm = _normalize(unit)
    if not norm or norm == "dimensionless":
        return float(value)
    if norm in _VOLUME_FACTORS:
        return convert_volume(value, unit)
    if norm in _CONCENTRATION_FACTORS:
        return convert_concentration(value, unit)
    if norm in _FLOW_FACTORS:
        return flow_to_l_per_day(value, unit)
    if norm in _SPECIAL_RATE_FACTORS:
        return float(value) * _SPECIAL_RATE_FACTORS[norm]
    if norm in {"1/day", "1/minute", "1/min", "1/hour", "1/h", "1/second", "1/sec", "1/s"}:
        return rate_to_per_day(value, unit)
    if norm.startswith("1/(mol") or norm.startswith("1/(micromolar"):
        return convert_kon(value, unit)
    ratio_prefixes = {
        "cell/day": lambda v: float(v),
        "cells/day": lambda v: float(v),
        "molecule/day": lambda v: float(v),
        "molecules/day": lambda v: float(v),
    }
    if norm in ratio_prefixes:
        return ratio_prefixes[norm](value)

    for prefix in ("1/day/", "1/hour/", "1/minute/", "1/second/"):
        if norm.startswith(prefix):
            denom = norm[len(prefix) :]
            converted = rate_to_per_day(value, prefix.rstrip("/"))
            return converted / _convert_denominator(denom)
    mole_prefixes = ("mole/", "mol/", "millimole/", "mmol/", "micromole/", "umol/")
    if norm.startswith(mole_prefixes):
        numer, denom = norm.split("/", 1)
        amount = convert_amount(value, numer, None)
        if denom:
            amount /= _convert_denominator(denom)
        return amount
    if norm.startswith("molecule/") or norm.startswith("molecules/"):
        _, denom = norm.split("/", 1)
        return float(value) / _convert_denominator(denom)
    if norm.startswith("1/") and norm[2:] in _VOLUME_FACTORS:
        return float(value) / convert_volume(1.0, norm[2:])
    if norm.startswith("1/(second"):
        return rate_to_per_day(value, "1/second")
    if norm in {"cell/milliliter", "cells/milliliter"}:
        return float(value) * 1000.0
    if norm in {"cell/liter", "cells/liter"}:
        return float(value)
    return float(value)


_MOLE_FACTORS = {
    "mole": 1.0,
    "mol": 1.0,
    "millimole": 1e-3,
    "mmol": 1e-3,
    "micromole": 1e-6,
    "umol": 1e-6,
}

_MASS_TO_GRAM = {
    "gram": 1.0,
    "g": 1.0,
    "milligram": 1e-3,
    "mg": 1e-3,
    "microgram": 1e-6,
    "ug": 1e-6,
}


def convert_amount(value: float, unit: str, molecular_weight_g_per_mol: float | None = None) -> float:
    norm = _normalize(unit)
    if not norm:
        return float(value)
    if norm in _MOLE_FACTORS:
        return float(value) * _MOLE_FACTORS[norm]
    if norm in _MASS_TO_GRAM:
        if molecular_weight_g_per_mol is None or molecular_weight_g_per_mol == 0.0:
            raise ValueError("Molecular weight (g/mol) required to convert mass to moles")
        grams = float(value) * _MASS_TO_GRAM[norm]
        return grams / molecular_weight_g_per_mol
    return float(value)


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
        return float(value) * factor
    return float(value)


_AREA_FACTORS = {
    "meter^2": 1.0,
    "metre^2": 1.0,
    "centimeter^2": 1e-4,
    "centimetre^2": 1e-4,
    "millimeter^2": 1e-6,
    "millimetre^2": 1e-6,
    "micrometer^2": 1e-12,
    "micrometre^2": 1e-12,
    "nanometer^2": 1e-18,
    "nanometre^2": 1e-18,
}


def convert_area(value: float, unit: str) -> float:
    norm = _normalize(unit)
    factor = _AREA_FACTORS.get(norm)
    if factor is not None:
        return float(value) * factor
    return float(value)


def _looks_like_concentration(units: str, dimension: Optional[str]) -> bool:
    dim_norm = _normalize(dimension) if dimension else ""
    if dim_norm:
        if any(token in dim_norm for token in ("concentration", "amount/vol", "mass/vol", "mol/vol")):
            return True
        if dim_norm in {"amount", "substance", "mass"}:
            return False
    unit_norm = _normalize(units)
    if not unit_norm:
        return False
    if "molar" in unit_norm or "mol/" in unit_norm or "mole/" in unit_norm:
        return True
    if "/l" in unit_norm or "/liter" in unit_norm or "/litre" in unit_norm:
        return True
    if "/" in unit_norm:
        return True
    return False


def normalise_dose_to_species(
    amount: float,
    amount_unit: str,
    *,
    species_units: Optional[str],
    species_dimension: Optional[str],
    compartment_volume_l: Optional[float],
    molecular_weight_g_per_mol: Optional[float] = None,
) -> Tuple[float, float]:
    """Convert an external dose into the target species units.

    Returns a tuple of (delta_in_species_units, amount_mol).
    """

    amount_mol = convert_amount(amount, amount_unit or "mole", molecular_weight_g_per_mol)
    if _looks_like_concentration(species_units or "", species_dimension):
        if compartment_volume_l in (None, 0.0):
            raise ValueError("Missing compartment volume for concentration species dose target")
        delta = amount_mol / float(compartment_volume_l)
        return delta, amount_mol
    return amount_mol, amount_mol


__all__ = [
    "convert_amount",
    "convert_area",
    "convert_concentration",
    "convert_kon",
    "convert_length",
    "convert_parameter_value",
    "convert_rate",
    "convert_volume",
    "flow_to_l_per_day",
    "normalise_dose_to_species",
    "rate_to_per_day",
    "time_to_day",
]  # pragma: no cover - used via import *
