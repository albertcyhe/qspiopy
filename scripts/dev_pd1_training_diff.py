#!/usr/bin/env python

"""Compute PD-1 white-box vs MATLAB diffs for selected training scenarios."""

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


def evaluate_scenario(
    scenario: str,
    frame: pd.DataFrame,
    params_path: Path,
    output_dir: Path,
) -> Tuple[float, float]:
    subset = frame[frame["scenario_id"] == scenario].sort_values("time_days")
    if subset.empty:
        raise ValueError(f"Scenario '{scenario}' not found in training parquet.")

    params = load_pd1_parameters_from_file(params_path)
    solver_cfg = SolverConfig(
        method="BDF",
        rtol=params.solver_rtol,
        atol=params.solver_atol,
        max_step=params.max_step_days,
        seed=None,
    )
    context = {
        "syn_pd1_pdl1": float(subset["syn_pd1_pdl1"].iloc[0]),
        "syn_pd1_pdl2": float(subset["syn_pd1_pdl2"].iloc[0]),
        "syn_pd1_apd1": float(subset["syn_pd1_apd1"].iloc[0]),
        "syn_pd1_apd1_pd1": float(subset["syn_pd1_apd1_pd1"].iloc[0]),
        "time_days": float(subset["time_days"].iloc[0]),
    }
    model = PD1WhiteboxModel.from_context(params, context, solver_config=solver_cfg)

    prev_time = float(subset["time_days"].iloc[0])
    preds_h: List[float] = []
    preds_block: List[float] = []
    rows = []
    first_index = subset.index[0]
    for idx, row in subset.iterrows():
        current_time = float(row["time_days"])
        dt = 0.0 if idx == first_index else max(current_time - prev_time, 0.0)
        outputs = model.step(float(row["drug_tumor_molar"]), dt)
        matlab_block = (float(row["syn_pd1_apd1"]) + 2.0 * float(row["syn_pd1_apd1_pd1"])) / max(
            model.total_pd1_density, 1e-12
        )
        preds_h.append(outputs.occupancy)
        preds_block.append(outputs.blocked_fraction)
        rows.append(
            {
                "time_days": current_time,
                "H_PD1_matlab": float(row["pd1_inhibition"]),
                "H_PD1_python": outputs.occupancy,
                "block_matlab": matlab_block,
                "block_python": outputs.blocked_fraction,
                "H_diff": outputs.occupancy - float(row["pd1_inhibition"]),
                "block_diff": outputs.blocked_fraction - matlab_block,
            }
        )
        prev_time = current_time

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"pd1_compare_{scenario}.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)

    matlab_h = subset["pd1_inhibition"].to_numpy()
    block_matlab = pd.DataFrame(rows)["block_matlab"].to_numpy()
    h_rmse = float(np.sqrt(np.mean((np.array(preds_h) - matlab_h) ** 2)))
    block_rmse = float(np.sqrt(np.mean((np.array(preds_block) - block_matlab) ** 2)))
    print(f"[pd1_diff] {scenario}: H_RMSE={h_rmse:.6g} block_RMSE={block_rmse:.6g} -> {out_path}")
    return h_rmse, block_rmse


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Diff PD-1 white-box against MATLAB training scenarios.")
    parser.add_argument("training_path", type=Path, help="Path to pd1_whitebox_training.parquet")
    parser.add_argument(
        "--parameter-file",
        type=Path,
        default=Path("parameters/example1_parameters.json"),
        help="Path to parameter JSON (default: parameters/example1_parameters.json)",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        default=["pd1_train_0004", "pd1_train_0582"],
        help="Scenario IDs to evaluate",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/dev"),
        help="Directory to write comparison CSVs",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    frame = pd.read_parquet(args.training_path)
    summary = []
    for scenario in args.scenarios:
        h_rmse, b_rmse = evaluate_scenario(scenario, frame, args.parameter_file, args.output_dir)
        summary.append((scenario, h_rmse, b_rmse))

    print("[pd1_diff] summary:")
    for scenario, h_rmse, b_rmse in summary:
        print(f"  {scenario}: H_RMSE={h_rmse:.6g} block_RMSE={b_rmse:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
