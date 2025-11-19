#!/usr/bin/env python

"""Reintegrate PD-1 white-box states to produce a clean ODE training set."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import pandas as pd
import numpy as np
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.offline.modules.pd1_params import load_pd1_parameters_from_file
from src.offline.modules.pd1_whitebox import PD1WhiteboxModel
from src.offline.segment_integrator import SolverConfig


def reintegrate_scenario(
    subset: pd.DataFrame,
    params_path: Path,
) -> pd.DataFrame:
    params = load_pd1_parameters_from_file(params_path)
    solver_cfg = SolverConfig(
        method="BDF",
        rtol=params.solver_rtol,
        atol=params.solver_atol,
        max_step=params.max_step_days,
        seed=None,
    )
    subset = subset.sort_values("time_days")
    first = subset.iloc[0]
    context = {
        "syn_pd1_pdl1": float(first.get("syn_pd1_pdl1", 0.0)),
        "syn_pd1_pdl2": float(first.get("syn_pd1_pdl2", 0.0)),
        "syn_pd1_apd1": float(first.get("syn_pd1_apd1", 0.0)),
        "syn_pd1_apd1_pd1": float(first.get("syn_pd1_apd1_pd1", 0.0)),
        "time_days": float(first["time_days"]),
    }
    model = PD1WhiteboxModel.from_context(params, context, solver_config=solver_cfg)

    rows: List[dict] = []
    prev_time = float(first["time_days"])
    for _, row in subset.iterrows():
        time = float(row["time_days"])
        dt = max(time - prev_time, 0.0)
        outputs = model.step(float(row["drug_tumor_molar"]), dt)
        rows.append(
            {
                "scenario_id": row["scenario_id"],
                "time_days": time,
                "drug_tumor_molar": row["drug_tumor_molar"],
                "drug_plasma_molar": row.get("drug_plasma_molar", np.nan),
                "H_PD1_python": outputs.occupancy,
                "block_fraction_python": outputs.blocked_fraction,
                "syn_pd1_pdl1_python": model.syn_pd1_pdl1,
                "syn_pd1_pdl2_python": model.syn_pd1_pdl2,
                "syn_pd1_apd1_python": model.syn_pd1_ab,
                "syn_pd1_apd1_pd1_python": model.syn_pd1_ab_pd1,
            }
        )
        prev_time = time

    return pd.DataFrame(rows)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Export clean PD-1 white-box ODE trajectories.")
    parser.add_argument(
        "training_path",
        type=Path,
        help="Existing training parquet providing scenarios and tumour drug traces.",
    )
    parser.add_argument(
        "--parameter-file",
        type=Path,
        default=Path("parameters/example1_parameters.json"),
        help="Path to PD-1 parameter JSON.",
    )
    parser.add_argument(
        "--output-parquet",
        type=Path,
        default=Path("artifacts/training/pd1_whitebox_training_clean.parquet"),
        help="Destination parquet for clean trajectories.",
    )
    parser.add_argument(
        "--max-scenarios",
        type=int,
        default=None,
        help="Limit number of scenarios processed (for quick runs).",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    frame = pd.read_parquet(args.training_path)
    scenario_ids = frame["scenario_id"].unique()
    if args.max_scenarios:
        scenario_ids = scenario_ids[: args.max_scenarios]

    outputs = []
    for scenario in scenario_ids:
        subset = frame[frame["scenario_id"] == scenario]
        if subset.empty:
            continue
        clean = reintegrate_scenario(subset, args.parameter_file)
        outputs.append(clean)
        print(f"[export_pd1_clean] processed {scenario} ({len(subset)} samples)")

    combined = pd.concat(outputs, ignore_index=True) if outputs else pd.DataFrame()
    out_path = args.output_parquet
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(out_path, index=False)
    print(f"[export_pd1_clean] wrote {out_path} ({len(combined)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
