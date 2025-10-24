"""Schema validation for frozen SimBiology snapshots."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
from pydantic import BaseModel, Field, ValidationError, field_validator

SNAPSHOT_FILES = [
    "configset.json",
    "species.csv",
    "parameters.csv",
    "compartments.csv",
    "rules.csv",
    "events.csv",
    "doses.csv",
    "variants.csv",
    "equations.txt",
]


def _clean_scalar(value):
    if pd.isna(value):
        return None
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "false"}:
            return lowered == "true"
        if lowered == "":
            return None
    return value


class SpeciesRow(BaseModel):
    identifier: str
    name: str
    compartment: str
    initial_amount: float
    initial_units: Optional[str] = None
    constant_amount: bool
    current_value: float
    boundary_condition: bool
    nonnegative: bool
    units: Optional[str] = None
    interpreted_dimension: str = Field(alias="interpreted_dimension")

    @field_validator("interpreted_dimension")
    def dimension_allowed(cls, value: str) -> str:
        allowed = {"amount", "concentration", "amount/vol", "amount/volume", "other", ""}
        if value not in allowed:
            raise ValueError(f"interpreted_dimension '{value}' is not recognised")
        return value


class ParameterRow(BaseModel):
    name: str
    value: float
    units: Optional[str] = None
    constant_value: bool
    initial_value: Optional[float] = None


class CompartmentRow(BaseModel):
    name: str
    capacity: float
    units: Optional[str] = None
    constant: bool


class RuleRow(BaseModel):
    rule_index: int
    name: Optional[str] = None
    type: str
    target: Optional[str] = None
    expression: Optional[str] = None

    @field_validator("type")
    def rule_type(cls, value: str) -> str:
        allowed = {"initialAssignment", "repeatedAssignment", "rate", "algebraic", "rateRule", "algebraicRule"}
        if value not in allowed:
            raise ValueError(f"unsupported rule type '{value}'")
        return value


class EventRow(BaseModel):
    event_index: int
    name: Optional[str] = None
    trigger: str
    delay: float = 0.0
    delay_type: Optional[str] = None
    assignments: Optional[str] = None


class DoseRow(BaseModel):
    dose_index: int
    name: str
    type: str
    target: str
    amount: float
    amount_units: Optional[str] = None
    start_time: float
    interval: float
    repeat_count: int


class VariantRow(BaseModel):
    variant_index: int
    name: str
    active: bool
    content: Optional[str] = None


def _validate_frame(frame: pd.DataFrame, model: BaseModel, filename: str) -> None:
    records = frame.to_dict(orient="records")
    for idx, record in enumerate(records, start=1):
        cleaned = {key: _clean_scalar(value) for key, value in record.items()}
        try:
            model.model_validate(cleaned)
        except ValidationError as exc:
            raise SystemExit(f"[schema] {filename} row {idx} invalid:\n{exc}") from exc


def _ensure_finite(frame: pd.DataFrame, columns: Iterable[str], filename: str) -> None:
    for column in columns:
        if column not in frame.columns:
            continue
        series = frame[column]
        for idx, value in enumerate(series, start=1):
            if pd.isna(value):
                continue
            if isinstance(value, (float, int)) and not math.isfinite(float(value)):
                raise SystemExit(f"[schema] {filename} row {idx} column '{column}' contains non-finite value {value!r}")


def _validate_config(config_path: Path) -> None:
    try:
        config = json.loads(config_path.read_text(encoding="utf8"))
    except json.JSONDecodeError as exc:  # pragma: no cover
        raise SystemExit(f"[schema] {config_path} is not valid JSON: {exc}") from exc
    required_keys = {
        "SolverType",
        "StopTime",
        "TimeUnits",
        "SolverOptions",
        "CompileOptions",
    }
    missing = required_keys.difference(config)
    if missing:
        raise SystemExit(f"[schema] configset.json missing keys: {sorted(missing)}")
    solver_opts = config.get("SolverOptions") or {}
    for key in ("AbsoluteTolerance", "RelativeTolerance"):
        if key not in solver_opts:
            raise SystemExit(f"[schema] SolverOptions missing '{key}'")
    compile_opts = config.get("CompileOptions") or {}
    if "DefaultSpeciesDimension" not in compile_opts:
        raise SystemExit("[schema] CompileOptions missing 'DefaultSpeciesDimension'")


def validate_snapshot(path: Path) -> None:
    path = path.resolve()
    if not path.is_dir():
        raise SystemExit(f"[schema] snapshot directory {path} is missing")

    missing = [name for name in SNAPSHOT_FILES if not (path / name).exists()]
    if missing:
        raise SystemExit(f"[schema] snapshot {path} missing files: {missing}")

    _validate_config(path / "configset.json")

    species = pd.read_csv(path / "species.csv")
    _validate_frame(species, SpeciesRow, "species.csv")
    _ensure_finite(species, ["initial_amount", "current_value"], "species.csv")

    parameters = pd.read_csv(path / "parameters.csv")
    _validate_frame(parameters, ParameterRow, "parameters.csv")

    compartments = pd.read_csv(path / "compartments.csv")
    _validate_frame(compartments, CompartmentRow, "compartments.csv")

    rules = pd.read_csv(path / "rules.csv")
    if not rules.empty:
        _validate_frame(rules, RuleRow, "rules.csv")

    events = pd.read_csv(path / "events.csv")
    if not events.empty:
        _validate_frame(events, EventRow, "events.csv")

    doses_path = path / "doses.csv"
    if doses_path.exists():
        doses = pd.read_csv(doses_path)
        if not doses.empty:
            _validate_frame(doses, DoseRow, "doses.csv")

    variants_path = path / "variants.csv"
    if variants_path.exists():
        variants = pd.read_csv(variants_path)
        if not variants.empty:
            _validate_frame(variants, VariantRow, "variants.csv")


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Validate schema for a frozen SimBiology snapshot directory")
    parser.add_argument("snapshot", type=Path, help="Path to snapshot directory (e.g. artifacts/matlab_frozen_model/example1)")
    args = parser.parse_args(argv)
    validate_snapshot(args.snapshot)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry
    raise SystemExit(main())
