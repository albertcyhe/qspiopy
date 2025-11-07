"""Generate unit-normalisation audit tables for snapshot parameters and doses."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Any

import pandas as pd

from src.offline.units import convert_parameter_value


def _load_rows(path: Path) -> list[dict[str, Any]]:
    if path.suffix.lower() == ".json":
        with path.open("r", encoding="utf8") as handle:
            data = json.load(handle)
        if isinstance(data, list):
            return data
        raise ValueError(f"Unsupported JSON structure in {path}")
    frame = pd.read_csv(path)
    return [row._asdict() for row in frame.itertuples(index=False)]


def audit_parameters(params_path: Path, output_path: Path) -> None:
    rows = _load_rows(params_path)
    records: list[dict[str, Any]] = []
    for row in rows:
        name = str(row.get("name"))
        unit = str(row.get("units", "") or "")
        value = row.get("value")
        if value is None:
            records.append({
                "name": name,
                "raw_value": value,
                "raw_unit": unit,
                "canonical_value": None,
                "note": "derived_only",
            })
            continue
        try:
            canonical = convert_parameter_value(float(value), unit)
            note = ""
        except Exception as exc:  # pragma: no cover - audit script only
            canonical = None
            note = f"ERROR: {exc}"
        records.append(
            {
                "name": name,
                "raw_value": value,
                "raw_unit": unit,
                "canonical_value": canonical,
                "note": note,
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(output_path, index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Audit unit conversions for snapshot parameters")
    parser.add_argument("parameters", type=Path, help="Path to parameters JSON/CSV exported by SimBiology")
    parser.add_argument("output", type=Path, help="Destination CSV for the audit table")
    args = parser.parse_args()
    audit_parameters(args.parameters, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
