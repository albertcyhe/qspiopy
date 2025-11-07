"""Emit a canonical units table from a frozen snapshot provenance map."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Any

import pandas as pd

from src.offline.snapshot import load_frozen_model


def _load_unit_map(provenance: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    raw = provenance.get("unit_normalisation_map") or provenance.get("parameter_graph")
    if not raw:
        return {}
    if isinstance(raw, str):
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            return {}
    if isinstance(raw, dict):
        return raw
    return {}


def emit_table(snapshot: str, output: Path) -> None:
    model = load_frozen_model(snapshot)
    unit_map = _load_unit_map(model.provenance)
    if not unit_map:
        raise RuntimeError(f"Snapshot '{snapshot}' does not expose a unit_normalisation_map")
    rows = []
    for name, meta in sorted(unit_map.items()):
        rows.append(
            {
                "name": name,
                "type": meta.get("type", ""),
                "unit": meta.get("unit", ""),
                "raw_value": meta.get("raw_value"),
                "canonical_value": meta.get("canonical_value"),
                "expression": meta.get("expression", ""),
                "dependencies": ";".join(meta.get("dependencies", [])),
                "note": meta.get("note", ""),
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Print the unit normalisation table for a snapshot")
    parser.add_argument("snapshot", help="Snapshot name or path (as accepted by load_frozen_model)")
    parser.add_argument("output", type=Path, help="Destination CSV path")
    args = parser.parse_args()
    emit_table(args.snapshot, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
