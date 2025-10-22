"""Compare the lightweight tutorial surrogate against the reference model."""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

from src.offline import simulate_reference, simulate_tutorial


@dataclass
class Scenario:
    name: str
    parameters: Tuple[Path, ...]
    therapy: str


def _compute_metrics(
    surrogate: pd.DataFrame, reference: pd.DataFrame, observable: str
) -> Dict[str, float]:
    merged = surrogate.merge(reference, on="time_days", suffixes=("_sur", "_ref"))

    y_true = merged[f"{observable}_ref"].to_numpy()
    y_pred = merged[f"{observable}_sur"].to_numpy()

    diff = y_pred - y_true
    rmse = float(np.sqrt(np.mean(diff**2)))
    ss_res = float(np.sum(diff**2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")

    auc_true = float(np.trapezoid(y_true, merged["time_days"]))
    auc_pred = float(np.trapezoid(y_pred, merged["time_days"]))
    delta_auc = float(auc_pred - auc_true)

    mean_abs = float(np.mean(np.abs(y_true))) or 1e-12
    rel_rmse = rmse / mean_abs
    max_pct = float(np.max(np.abs(diff) / (np.abs(y_true) + 1e-12)))

    return {
        "rmse": rmse,
        "r2": r2,
        "delta_auc": delta_auc,
        "relative_rmse": rel_rmse,
        "max_fractional_error": max_pct,
    }


def _tgi(volume_control: pd.DataFrame, volume_treated: pd.DataFrame, column: str) -> float:
    ctrl = volume_control.iloc[-1][column]
    trt = volume_treated.iloc[-1][column]
    if ctrl <= 0:
        return 0.0
    return 1.0 - float(trt / ctrl)


def _run_scenario(scenario: Scenario, output_dir: Path) -> Dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    surrogate = simulate_tutorial(scenario.parameters, therapy=scenario.therapy)
    reference = simulate_reference(scenario.parameters, therapy=scenario.therapy)

    surrogate_path = output_dir / f"{scenario.name}_surrogate.csv"
    reference_path = output_dir / f"{scenario.name}_reference.csv"

    surrogate.to_frame().to_csv(surrogate_path, index=False)
    reference.to_frame().to_csv(reference_path, index=False)

    return {"surrogate": surrogate_path, "reference": reference_path}


def _performance_benchmark(parameters: Tuple[Path, ...], replicates: int) -> Dict[str, float]:
    start = time.perf_counter()
    for _ in range(replicates):
        simulate_reference(parameters, therapy="anti_pd1")
    ref_time = time.perf_counter() - start

    start = time.perf_counter()
    for _ in range(replicates):
        simulate_tutorial(parameters, therapy="anti_pd1")
    sur_time = time.perf_counter() - start

    return {"reference_seconds": ref_time, "surrogate_seconds": sur_time}


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate the Python surrogate against the reference equations")
    parser.add_argument("--output", type=Path, default=Path("artifacts/validation"))
    parser.add_argument("--replicates", type=int, default=50, help="Number of runs for the performance benchmark")
    args = parser.parse_args(argv)

    scenarios = [
        Scenario("example1_control", (Path("parameters/example1_parameters.json"),), "none"),
        Scenario("example1_treated", (Path("parameters/example1_parameters.json"),), "anti_pd1"),
        Scenario("example2_treated", (Path("parameters/example2_parameters.json"),), "anti_pd1"),
    ]

    results: List[Dict[str, Path]] = []
    metrics_records: List[Dict[str, object]] = []

    for scenario in scenarios:
        paths = _run_scenario(scenario, args.output)
        results.append({"scenario": scenario.name, **paths})

        surrogate = pd.read_csv(paths["surrogate"])
        reference = pd.read_csv(paths["reference"])

        for observable in ["tumour_volume_l", "pd1_occupancy", "tcell_density_per_ul"]:
            metrics = _compute_metrics(surrogate, reference, observable)
            metrics_records.append(
                {
                    "scenario": scenario.name,
                    "observable": observable,
                    **metrics,
                }
            )

    example1_sur = pd.read_csv(args.output / "example1_treated_surrogate.csv")
    example1_ref = pd.read_csv(args.output / "example1_treated_reference.csv")
    example1_sur_ctrl = pd.read_csv(args.output / "example1_control_surrogate.csv")
    example1_ref_ctrl = pd.read_csv(args.output / "example1_control_reference.csv")

    tgi_sur = _tgi(example1_sur_ctrl, example1_sur, "tumour_volume_l")
    tgi_ref = _tgi(example1_ref_ctrl, example1_ref, "tumour_volume_l")
    metrics_records.append(
        {
            "scenario": "example1_treated",
            "observable": "tgi",
            "rmse": abs(tgi_sur - tgi_ref),
            "r2": float("nan"),
            "delta_auc": tgi_sur - tgi_ref,
        }
    )

    metrics_df = pd.DataFrame(metrics_records)
    metrics_path = args.output / "metrics.csv"
    metrics_df.to_csv(metrics_path, index=False)

    benchmark = _performance_benchmark((Path("parameters/example1_parameters.json"),), args.replicates)
    with (args.output / "performance.json").open("w", encoding="utf8") as handle:
        json.dump({"replicates": args.replicates, **benchmark}, handle, indent=2)

    registry_df = pd.DataFrame(results)
    registry_df.to_csv(args.output / "artefacts.csv", index=False)

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry
    raise SystemExit(main())

