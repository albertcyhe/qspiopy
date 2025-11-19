#!/usr/bin/env python

"""Compare MATLAB probe CSVs against the PD-1 white-box integration."""

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

from src.offline.modules.pd1_params import load_pd1_parameters_from_file
from src.offline.modules.pd1_whitebox import PD1WhiteboxModel
from src.offline.segment_integrator import SolverConfig


def evaluate_probe(csv_path: Path, params_path: Path) -> Tuple[float, float, float]:
    frame = pd.read_csv(csv_path)
    params = load_pd1_parameters_from_file(params_path)
    solver_cfg = SolverConfig(
        method="BDF",
        rtol=params.solver_rtol,
        atol=params.solver_atol,
        max_step=params.max_step_days,
        seed=None,
    )

    first = frame.iloc[0]
    context = {
        "syn_pd1_pdl1": float(first["PD1_PDL1"]),
        "syn_pd1_pdl2": float(first["PD1_PDL2"]),
        "syn_pd1_apd1": float(first["PD1_aPD1"]),
        "syn_pd1_apd1_pd1": float(first["PD1_aPD1_PD1"]),
        "time_days": float(first["time_days"]),
    }
    model = PD1WhiteboxModel.from_context(params, context, solver_config=solver_cfg)

    prev_time = float(first["time_days"])
    sims_pdl1: List[float] = []
    sims_ab: List[float] = []
    sims_dimer: List[float] = []
    blocks: List[float] = []
    for _, row in frame.iterrows():
        time = float(row["time_days"])
        dt = max(time - prev_time, 0.0)
        model.step(float(row["aPD1_molar"]), dt)
        sims_pdl1.append(model.syn_pd1_pdl1)
        sims_ab.append(model.syn_pd1_ab)
        sims_dimer.append(model.syn_pd1_ab_pd1)
        blocks.append(model._blocked_fraction())
        prev_time = time

    matlab_block = (frame["PD1_aPD1"] + 2.0 * frame["PD1_aPD1_PD1"]) / max(params.pd1_surface_density(), 1e-12)

    def rmse(pred: np.ndarray, target: np.ndarray) -> float:
        return float(np.sqrt(np.mean((pred - target) ** 2)))

    rmse_pdl1 = rmse(np.array(sims_pdl1), frame["PD1_PDL1"].to_numpy())
    rmse_block = rmse(np.array(blocks), matlab_block.to_numpy())
    rmse_ab = rmse(np.array(sims_ab) + 2.0 * np.array(sims_dimer), (frame["PD1_aPD1"] + 2.0 * frame["PD1_aPD1_PD1"]).to_numpy())
    return rmse_pdl1, rmse_ab, rmse_block


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="PD-1 probe diff against white-box integration.")
    parser.add_argument(
        "probe_csvs",
        nargs="+",
        help="Probe CSV paths (from matlab/scripts/dev_pd1_training_probe.m).",
    )
    parser.add_argument(
        "--parameter-file",
        type=Path,
        default=Path("parameters/example1_parameters.json"),
        help="PD-1 parameter JSON path.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    for path_str in args.probe_csvs:
        csv_path = Path(path_str)
        rmse_pdl1, rmse_bound, rmse_block = evaluate_probe(csv_path, args.parameter_file)
        print(
            f"[probe_diff] {csv_path.name}: "
            f"RMSE_PD1_PDL1={rmse_pdl1:.6g}, "
            f"RMSE_bound={rmse_bound:.6g}, "
            f"RMSE_block={rmse_block:.6g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
