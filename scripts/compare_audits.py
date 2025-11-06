"""Dose and event audit comparison utilities."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd


def _canonicalise_dose(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    rename_map = {
        "delta_applied": "delta_state_value",
        "delta_value": "delta_state_value",
        "amount_mg": "mg",
        "compartment": "compartment_name",
        "Vc_L": "compartment_volume_l",
        "dose_name": "name",
        "delta_amount": "delta_amount_mol",
    }
    df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns}, inplace=True)
    if "time_days" in df.columns and "time_hours" not in df.columns:
        df["time_hours"] = df["time_days"] * 24.0
    if "compartment_volume_l" not in df.columns:
        df["compartment_volume_l"] = np.nan
    if "mg" not in df.columns and "amount_mg" in df.columns:
        df["mg"] = df["amount_mg"]
    if "delta_state_value" not in df.columns:
        if "delta_applied" in df.columns:
            df["delta_state_value"] = df["delta_applied"]
        else:
            df["delta_state_value"] = np.nan
    if "delta_amount_mol" not in df.columns:
        if "amount_moles" in df.columns:
            df["delta_amount_mol"] = df["amount_moles"]
        else:
            df["delta_amount_mol"] = np.nan
    if "amount_moles" not in df.columns and "delta_amount_mol" in df.columns:
        df["amount_moles"] = df["delta_amount_mol"]
    if "interpreted_dimension" not in df.columns:
        df["interpreted_dimension"] = ""
    if "units" not in df.columns:
        if "species_units" in df.columns:
            df["units"] = df["species_units"]
        else:
            df["units"] = ""
    df["dimension_key"] = df.apply(
        lambda row: row["units"] if str(row["units"]) else str(row["interpreted_dimension"]), axis=1
    )
    df = df.sort_values(["time_hours", "target"]).reset_index(drop=True)
    return df


def _canonicalise_events(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "time_fire" in df.columns and "time_hours" not in df.columns:
        df["time_hours"] = df["time_fire"] * 24.0
    if "index_in_model" not in df.columns and "event_index" in df.columns:
        df["index_in_model"] = df["event_index"]
    df = df.sort_values(["time_hours", "index_in_model"]).reset_index(drop=True)
    return df


def _near(a: float, b: float, rtol: float = 1e-6, atol: float = 1e-12) -> bool:
    return abs(a - b) <= atol + rtol * abs(b)


def compare_doses(py_path: Path, ref_path: Path) -> List[str]:
    py = _canonicalise_dose(pd.read_csv(py_path))
    rf = _canonicalise_dose(pd.read_csv(ref_path))
    issues: List[str] = []
    if len(py) != len(rf):
        issues.append(f"dose count mismatch: py={len(py)} ref={len(rf)}")
    n = min(len(py), len(rf))
    for idx in range(n):
        row_p = py.iloc[idx]
        row_r = rf.iloc[idx]
        if not _near(row_p["time_hours"], row_r["time_hours"]) or row_p["target"] != row_r["target"]:
            issues.append(
                f"[{idx}] time/target mismatch: "
                f"{row_p['time_hours']}/{row_p['target']} vs {row_r['time_hours']}/{row_r['target']}"
            )
            continue
        dimension_p = row_p["dimension_key"]
        dimension_r = row_r["dimension_key"]
        if dimension_p != dimension_r:
            issues.append(
                f"[{idx}] dimension mismatch: {dimension_p} vs {dimension_r}"
            )
        if not _near(row_p["delta_state_value"], row_r["delta_state_value"]):
            issues.append(
                f"[{idx}] delta_state mismatch: py={row_p['delta_state_value']:.6g} "
                f"ref={row_r['delta_state_value']:.6g}"
            )
        if not _near(row_p["mg"], row_r["mg"]):
            issues.append(f"[{idx}] mg mismatch: py={row_p['mg']:.6g} ref={row_r['mg']:.6g}")
        delta_amount_p = row_p.get("delta_amount_mol", np.nan)
        delta_amount_r = row_r.get("delta_amount_mol", np.nan)
        if not (np.isnan(delta_amount_p) or np.isnan(delta_amount_r)):
            if not _near(delta_amount_p, delta_amount_r):
                issues.append(
                    f"[{idx}] delta_amount mismatch: py={delta_amount_p:.6g} ref={delta_amount_r:.6g}"
                )
        compartment_vc = row_p.get("compartment_volume_l", np.nan)
        if isinstance(compartment_vc, (float, int)) and not np.isnan(compartment_vc):
            dim_lower = str(row_p["interpreted_dimension"]).lower()
            if "conc" in dim_lower or "/l" in dim_lower or "mol" in dim_lower:
                state_to_amount_p = row_p["delta_state_value"] * row_p["compartment_volume_l"]
                state_to_amount_r = row_r["delta_state_value"] * row_r["compartment_volume_l"]
                if not _near(row_p["delta_amount_mol"], state_to_amount_p):
                    issues.append(
                        f"[{idx}] surrogate Δamount mismatch: {row_p['delta_amount_mol']:.6g} vs "
                        f"{state_to_amount_p:.6g} (=Δstate*Vc)"
                    )
                if not _near(row_r["delta_amount_mol"], state_to_amount_r):
                    issues.append(
                        f"[{idx}] reference Δamount mismatch: {row_r['delta_amount_mol']:.6g} vs "
                        f"{state_to_amount_r:.6g} (=Δstate*Vc)"
                    )
            else:
                if not np.isnan(row_p["delta_amount_mol"]) and not _near(
                    row_p["delta_amount_mol"], row_p["delta_state_value"]
                ):
                    issues.append(
                        f"[{idx}] surrogate Δstate ({row_p['delta_state_value']:.6g}) "
                        f"!= Δamount ({row_p['delta_amount_mol']:.6g}) for mass target"
                    )
                if not np.isnan(row_r["delta_amount_mol"]) and not _near(
                    row_r["delta_amount_mol"], row_r["delta_state_value"]
                ):
                    issues.append(
                        f"[{idx}] reference Δstate ({row_r['delta_state_value']:.6g}) "
                        f"!= Δamount ({row_r['delta_amount_mol']:.6g}) for mass target"
                    )
    return issues


def compare_events(py_evt_path: Path, ref_evt_path: Path, tol_h: float = 0.5) -> List[str]:
    py_events = _canonicalise_events(pd.read_csv(py_evt_path))
    ref_events = _canonicalise_events(pd.read_csv(ref_evt_path))
    issues: List[str] = []
    if len(py_events) != len(ref_events):
        issues.append(f"event count mismatch: py={len(py_events)} ref={len(ref_events)}")
    n = min(len(py_events), len(ref_events))
    for idx in range(n):
        row_p = py_events.iloc[idx]
        row_r = ref_events.iloc[idx]
        dt = abs(float(row_p["time_hours"]) - float(row_r["time_hours"]))
        if dt > tol_h:
            issues.append(f"[{idx}] event time Δt={dt:.3f} h > {tol_h} h")
        if row_p["index_in_model"] != row_r["index_in_model"]:
            issues.append(
                f"[{idx}] event order mismatch: {row_p['index_in_model']} vs {row_r['index_in_model']}"
            )
    return issues


def main(argv: List[str]) -> int:
    if len(argv) != 4:
        print(
            "Usage: python -m scripts.compare_audits "
            "<surrogate_dose_audit.csv> <reference_dose_audit.csv> "
            "<surrogate_events.csv> <reference_events.csv>"
        )
        return 1

    py_dose, ref_dose, py_evt, ref_evt = map(Path, argv)
    issues = compare_doses(py_dose, ref_dose)
    issues += compare_events(py_evt, ref_evt)

    if issues:
        print("❌ Audit mismatches:")
        for msg in issues:
            print(" -", msg)
        return 2

    print("✅ Dose & event audits aligned.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
