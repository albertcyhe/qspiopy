#!/usr/bin/env python

"""Compare MATLAB T-cell probe outputs against the Python white-box evolution."""

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
from src.offline.modules.tcell_whitebox import TCellWhiteboxModel


def _rmse(pred: np.ndarray, target: np.ndarray) -> float:
    return float(np.sqrt(np.mean((pred - target) ** 2)))


def _seed_context_from_probe(frame: pd.DataFrame, model_name: str) -> Tuple[dict, dict]:
    first = frame.iloc[0]
    model = load_frozen_model(model_name)
    parameters = dict(model.parameters or {})
    compartments = model.compartments or {}
    volume_l = float(first["tumour_volume_l"])
    context = {
        "time_days": float(first["time_days"]),
        "tumour_volume_l": volume_l,
        "V_T": volume_l * 1e6,
        "V_T.T1": float(first["Tumor_T1"]),
        "V_T.T0": float(first["Tumor_T0"]),
        "V_T.T_exh": float(first["Tumor_Texh"]),
        "T_exh": float(first["Tumor_Texh"]),
        "tcell_density_per_ul": float(first["tcell_density_per_ul"]),
        "C1": float(first["Tumor_C1"]),
        "C_total": float(first["Tumor_C1"]),
        "V_T.C1": float(first["Tumor_C1"]),
        "H_PD1_C1": float(first["H_PD1_C1"]),
        "H_mAPC": float(first.get("H_mAPC", 1.0)),
        "H_P1": float(first.get("H_P1", 1.0)),
        "V_C": float(compartments.get("V_C", 5.0)),
        "V_P": float(compartments.get("V_P", 60.0)),
        "V_LN": float(compartments.get("V_LN", 1.0)),
        "cell": float(parameters.get("cell", 1.0)),
    }
    context["V_LN.nT1"] = float(first["LN_nT1"])
    context["V_LN.aT1"] = float(first["LN_aT1"])
    context["V_LN.T1"] = float(first["LN_T1"])
    context["V_C.nT1"] = float(first["C_nT1"])
    context["V_C.T1"] = float(first["C_T1"])
    context["V_P.nT1"] = float(first["P_nT1"])
    context["V_P.T1"] = float(first["P_T1"])
    return context, parameters


def _suggest_density_scale(frame: pd.DataFrame) -> float:
    mask = (frame["Tumor_T1"] > 0.0) & (frame["tumour_volume_l"] > 0.0)
    if not mask.any():
        return float("nan")
    subset = frame.loc[mask]
    ratio = subset["tcell_density_per_ul"] * subset["tumour_volume_l"] * 1e6 / subset["Tumor_T1"]
    ratio = ratio[np.isfinite(ratio)]
    if not len(ratio):
        return float("nan")
    return float(np.median(ratio))


def integrate_from_probe(
    csv_path: Path,
    snapshot_name: str,
    density_scale: float | None = None,
    *,
    track_derivs: bool = False,
) -> Tuple[float, float, float]:
    frame = pd.read_csv(csv_path)
    context, parameters = _seed_context_from_probe(frame, snapshot_name)
    if density_scale is not None and np.isfinite(density_scale):
        parameters["tcell_density_scale"] = float(density_scale)
    tcell_model = TCellWhiteboxModel.from_context(parameters, context)

    pred_density: List[float] = []
    pred_t_eff: List[float] = []
    pred_d_t_eff: List[float] = []

    prev_time = float(context["time_days"])
    for _, row in frame.iterrows():
        current_time = float(row["time_days"])
        dt = max(current_time - prev_time, 0.0)
        prev_time = current_time
        volume_l = float(row["tumour_volume_l"])
        context["tumour_volume_l"] = volume_l
        context["V_T"] = volume_l * 1e6
        context["time_days"] = current_time
        context["V_T.T1"] = float(row.get("Tumor_T1", context.get("V_T.T1", 0.0)))
        context["V_T.T0"] = float(row.get("Tumor_T0", context.get("V_T.T0", 0.0)))
        context["V_T.T_exh"] = float(row.get("Tumor_Texh", context.get("V_T.T_exh", 0.0)))
        context["T_exh"] = context["V_T.T_exh"]
        context["C1"] = float(row.get("Tumor_C1", context.get("C1", 0.0)))
        context["C_total"] = context["C1"]
        context["V_T.C1"] = context["C1"]
        context["V_T.T0"] = float(row.get("Tumor_T0", context.get("V_T.T0", 0.0)))
        context["V_T.T_exh"] = float(row.get("Tumor_Texh", context.get("V_T.T_exh", 0.0)))
        context["H_PD1_C1"] = float(row.get("H_PD1_C1", context.get("H_PD1_C1", 0.0)))
        context["H_mAPC"] = float(row.get("H_mAPC", context.get("H_mAPC", 1.0)))
        context["H_P1"] = float(row.get("H_P1", context.get("H_P1", 1.0)))
        occ = float(context.get("H_PD1_C1", 0.0))
        tcell_model.step(context, dt, occ)
        pred_density.append(context.get("tcell_density_per_ul", 0.0))
        pred_t_eff.append(context.get("V_T.T1", 0.0))
        rhs = getattr(tcell_model, "debug_last_rhs", None)
        if rhs is not None and len(rhs) > 7:
            pred_d_t_eff.append(float(rhs[7]))
        else:
            pred_d_t_eff.append(float("nan"))

    rmse_density = _rmse(np.array(pred_density), frame["tcell_density_per_ul"].to_numpy())
    rmse_eff = _rmse(np.array(pred_t_eff), frame["Tumor_T1"].to_numpy())
    if track_derivs:
        pred_arr = np.array(pred_d_t_eff)
        target_arr = frame["d_Tumor_T1_dt"].to_numpy()
        mask = np.isfinite(pred_arr) & np.isfinite(target_arr)
        rmse_d_t_eff = _rmse(pred_arr[mask], target_arr[mask]) if mask.any() else float("nan")
    else:
        rmse_d_t_eff = float("nan")
    return rmse_density, rmse_eff, rmse_d_t_eff


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Diff MATLAB T-cell probe against Python integration.")
    parser.add_argument("probe_csvs", nargs="+", help="Paths to tcell_probe CSVs.")
    parser.add_argument(
        "--snapshot-dir",
        type=Path,
        default=Path("artifacts/matlab_frozen_model/example1"),
        help="Snapshot directory containing the frozen snapshot.",
    )
    parser.add_argument(
        "--estimate-density-scale",
        action="store_true",
        help="Estimate the tcell_density_scale factor from the probe CSV.",
    )
    parser.add_argument(
        "--density-scale",
        type=float,
        default=None,
        help="Override tcell_density_scale when integrating the probe (useful to test calibration).",
    )
    parser.add_argument(
        "--track-derivs",
        action="store_true",
        help="Also report RMSE for d(Tumor_T1)/dt to debug RHS contributions.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    snapshot_name = str(args.snapshot_dir)
    for path_str in args.probe_csvs:
        csv_path = Path(path_str)
        rmse_density, rmse_eff, rmse_d_t_eff = integrate_from_probe(
            csv_path,
            snapshot_name,
            args.density_scale,
            track_derivs=args.track_derivs,
        )
        print(
            f"[tcell_probe_diff] {csv_path.name}: "
            f"RMSE_density={rmse_density:.6g}, "
            f"RMSE_tumour_T1={rmse_eff:.6g}"
            + (f", RMSE_dTumor_T1_dt={rmse_d_t_eff:.6g}" if args.track_derivs else "")
        )
        if args.estimate_density_scale:
            frame = pd.read_csv(csv_path)
            estimate = _suggest_density_scale(frame)
            if np.isfinite(estimate):
                print(f"  -> Suggested tcell_density_scale ≈ {estimate:.6g}")
            else:
                print("  -> Could not estimate tcell_density_scale (insufficient data).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
