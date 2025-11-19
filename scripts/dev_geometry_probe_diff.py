#!/usr/bin/env python

"""Compare MATLAB geometry probe outputs against the Python white-box module."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.offline import load_frozen_model
from src.offline.modules.geometry_whitebox import GeometryWhiteboxModel


def _rmse(pred: np.ndarray, target: np.ndarray) -> float:
    return float(np.sqrt(np.mean((pred - target) ** 2)))


def _seed_geometry_context(frame: pd.DataFrame, snapshot_name: str) -> Tuple[dict, dict]:
    first = frame.iloc[0]
    model = load_frozen_model(snapshot_name)
    parameters = model.parameters or {}
    volume_l = float(first["tumour_volume_l"])
    context = {
        "time_days": float(first["time_days"]),
        "tumour_volume_l": volume_l,
        "V_T": volume_l * 1e6,
        "V_T.C1": float(first["Tumor_C1"]),
        "V_T.C_x": float(first.get("Tumor_Cx", first["Tumor_C1"] * 0.01)),
        "C1": float(first["Tumor_C1"]),
        "C_x": float(first.get("Tumor_Cx", first["Tumor_C1"] * 0.01)),
        "V_T.T1": float(first["Tumor_T1"]),
        "V_T.T0": float(first["Tumor_T0"]),
        "V_T.T_exh": float(first["Tumor_Texh"]),
        "T_exh": float(first["Tumor_Texh"]),
        "H_PD1_C1": float(first["H_PD1_C1"]),
        "cell": float(parameters.get("cell", 1.0)),
    }
    return context, parameters


def integrate_geometry_from_probe(csv_path: Path, snapshot_name: str) -> Tuple[float, float, float]:
    frame = pd.read_csv(csv_path)
    context, parameters = _seed_geometry_context(frame, snapshot_name)
    geom_model = GeometryWhiteboxModel.from_context(parameters, context)

    pred_volume: List[float] = []
    pred_c1: List[float] = []
    pred_cx: List[float] = []

    prev_time = float(context["time_days"])
    first_row = True
    for _, row in frame.iterrows():
        current_time = float(row["time_days"])
        dt = max(current_time - prev_time, 0.0)
        prev_time = current_time
        context["time_days"] = current_time
        context["V_T.T1"] = float(row["Tumor_T1"])
        context["V_T.T0"] = float(row["Tumor_T0"])
        context["V_T.T_exh"] = float(row["Tumor_Texh"])
        context["T_exh"] = float(row["Tumor_Texh"])
        context["H_PD1_C1"] = float(row["H_PD1_C1"])
        if first_row:
            volume_l = float(row["tumour_volume_l"])
            context["tumour_volume_l"] = volume_l
            context["V_T"] = volume_l * 1e6
            context["V_T.C1"] = float(row["Tumor_C1"])
            context["C1"] = float(row["Tumor_C1"])
            context["V_T.C_x"] = float(row.get("Tumor_Cx", context.get("V_T.C_x", 0.0)))
            context["C_x"] = context["V_T.C_x"]
            first_row = False
        occ = float(context["H_PD1_C1"])
        geom_model.step(context, dt, occ)
        pred_volume.append(context.get("tumour_volume_l", 0.0))
        pred_c1.append(context.get("V_T.C1", 0.0))
        pred_cx.append(context.get("V_T.C_x", 0.0))

    rmse_volume = _rmse(np.array(pred_volume), frame["tumour_volume_l"].to_numpy())
    rmse_c1 = _rmse(np.array(pred_c1), frame["Tumor_C1"].to_numpy())
    if "Tumor_Cx" in frame.columns:
        rmse_cx = _rmse(np.array(pred_cx), frame["Tumor_Cx"].to_numpy())
    else:
        rmse_cx = float("nan")
    return rmse_volume, rmse_c1, rmse_cx


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Diff MATLAB geometry probe CSVs against the Python white-box.")
    parser.add_argument("probe_csvs", nargs="+", help="Paths to geometry probe CSVs.")
    parser.add_argument(
        "--snapshot-dir",
        type=Path,
        default=Path("artifacts/matlab_frozen_model/example1"),
        help="Snapshot directory containing the frozen snapshot.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    snapshot_name = str(args.snapshot_dir)
    for path_str in args.probe_csvs:
        csv_path = Path(path_str)
        rmse_volume, rmse_c1, rmse_cx = integrate_geometry_from_probe(csv_path, snapshot_name)
        msg = (
            f"[geometry_probe_diff] {csv_path.name}: "
            f"RMSE_volume={rmse_volume:.6g}, "
            f"RMSE_C1={rmse_c1:.6g}"
        )
        if np.isfinite(rmse_cx):
            msg += f", RMSE_Cx={rmse_cx:.6g}"
        print(msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
