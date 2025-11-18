#!/usr/bin/env python
"""Compare PD-1 synapse derivatives vs MATLAB finite differences."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.offline.modules.pd1_params import load_pd1_parameters_from_file
from src.offline.modules.pd1_whitebox import pd1_synapse_reaction_terms, pd1_synapse_rhs


STATE_COLUMNS = ["syn_pd1_pdl1", "syn_pd1_pdl2", "syn_pd1_apd1", "syn_pd1_apd1_pd1"]


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare PD-1 RHS vs MATLAB finite-difference derivatives.")
    parser.add_argument("training_path", type=Path, help="Path to pd1_whitebox_training.parquet")
    parser.add_argument(
        "--parameter-file",
        type=Path,
        default=Path("parameters/example1_parameters.json"),
        help="Parameter JSON (default: parameters/example1_parameters.json)",
    )
    parser.add_argument(
        "--scenarios",
        nargs="+",
        default=["pd1_train_0004"],
        help="Scenario IDs to evaluate (default: pd1_train_0004).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/dev"),
        help="Directory to write derivative comparison CSVs.",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def build_python_rhs(
    params,
    *,
    total_pd1: float,
    total_pdl1: float,
    total_pdl2: float,
    ab_molar: float,
    state_vector: np.ndarray,
) -> np.ndarray:
    ab_effective = max(float(ab_molar), 0.0) / max(params.gamma_c_nivolumab, 1e-9)
    return pd1_synapse_rhs(
        state_vector,
        params,
        total_pd1_density=total_pd1,
        total_pdl1_density=total_pdl1,
        total_pdl2_density=total_pdl2,
        ab_effective_molar=ab_effective,
    )


def finite_difference(values: np.ndarray, times: np.ndarray, idx: int) -> float:
    dt = times[idx + 1] - times[idx - 1]
    if dt == 0.0:
        return 0.0
    return (values[idx + 1] - values[idx - 1]) / dt


def evaluate_scenario(frame: pd.DataFrame, scenario: str, params, output_dir: Path) -> None:
    subset = frame[frame["scenario_id"] == scenario].sort_values("time_days").reset_index(drop=True)
    if len(subset) < 3:
        raise ValueError(f"Scenario '{scenario}' requires at least 3 time points for central diff.")

    total_pd1 = params.pd1_surface_density()
    total_pdl1 = params.pdl1_surface_density()
    total_pdl2 = params.pdl2_surface_density()

    times = subset["time_days"].to_numpy(dtype=float)
    state_matrix = subset[STATE_COLUMNS].to_numpy(dtype=float)
    drug_trace = subset["drug_tumor_molar"].to_numpy(dtype=float)

    records: List[dict] = []
    for idx in range(1, len(subset) - 1):
        time = times[idx]
        diffs_matlab = {}
        for col in range(len(STATE_COLUMNS)):
            diffs_matlab[STATE_COLUMNS[col]] = finite_difference(state_matrix[:, col], times, idx)

        rhs = build_python_rhs(
            params,
            total_pd1=total_pd1,
            total_pdl1=total_pdl1,
            total_pdl2=total_pdl2,
            ab_molar=float(drug_trace[idx]),
            state_vector=state_matrix[idx, :],
        )
        reaction_89, reaction_90, reaction_91, reaction_92, pd1_free, pdl1_free, pdl2_free = pd1_synapse_reaction_terms(
            state_matrix[idx, :],
            params,
            total_pd1_density=total_pd1,
            total_pdl1_density=total_pdl1,
            total_pdl2_density=total_pdl2,
            ab_effective_molar= max(float(drug_trace[idx]), 0.0) / max(params.gamma_c_nivolumab, 1e-9),
        )
        records.append(
            {
                "time_days": time,
                "drug_tumor_molar": drug_trace[idx],
                "dpdl1_python": rhs[0],
                "dpdl1_matlab": diffs_matlab["syn_pd1_pdl1"],
                "dpdl1_diff": rhs[0] - diffs_matlab["syn_pd1_pdl1"],
                "dpdl2_python": rhs[1],
                "dpdl2_matlab": diffs_matlab["syn_pd1_pdl2"],
                "dpdl2_diff": rhs[1] - diffs_matlab["syn_pd1_pdl2"],
                "dab_python": rhs[2],
                "dab_matlab": diffs_matlab["syn_pd1_apd1"],
                "dab_diff": rhs[2] - diffs_matlab["syn_pd1_apd1"],
                "ddimer_python": rhs[3],
                "ddimer_matlab": diffs_matlab["syn_pd1_apd1_pd1"],
                "ddimer_diff": rhs[3] - diffs_matlab["syn_pd1_apd1_pd1"],
                "reaction_89": reaction_89,
                "reaction_90": reaction_90,
                "reaction_91": reaction_91,
                "reaction_92": reaction_92,
                "pd1_free": pd1_free,
                "pdl1_free": pdl1_free,
                "pdl2_free": pdl2_free,
            }
        )

    df = pd.DataFrame(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"pd1_derivative_diff_{scenario}.csv"
    df.to_csv(out_path, index=False)
    max_abs = df[["dpdl1_diff", "dpdl2_diff", "dab_diff", "ddimer_diff"]].abs().max()
    print(f"[pd1_derivative_diff] {scenario}: max|diff| -> "
          f"dPdl1={max_abs['dpdl1_diff']:.3e}, "
          f"dPdl2={max_abs['dpdl2_diff']:.3e}, "
          f"dAb={max_abs['dab_diff']:.3e}, "
          f"dDimer={max_abs['ddimer_diff']:.3e} ({out_path})")


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    frame = pd.read_parquet(args.training_path)
    params = load_pd1_parameters_from_file(args.parameter_file)
    for scenario in args.scenarios:
        evaluate_scenario(frame, scenario, params, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
