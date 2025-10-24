"""Compare MATLAB SimBiology trajectories against the Python reference model."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Scenario:
    name: str
    has_control: bool
    observables: Sequence[str] = ("tumour_volume_l", "pd1_occupancy", "tcell_density_per_ul")


MATLAB_SUFFIX = "_matlab"
PYTHON_SUFFIX = "_python"


def _load_frame(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def _compute_metrics(matlab: pd.DataFrame, python: pd.DataFrame, observable: str) -> dict[str, float]:
    merged = matlab.merge(python, on="time_days", suffixes=(MATLAB_SUFFIX, PYTHON_SUFFIX))
    if merged.empty:
        raise ValueError(f"No overlapping time points for observable {observable}")

    true = merged[f"{observable}{MATLAB_SUFFIX}"].to_numpy()
    pred = merged[f"{observable}{PYTHON_SUFFIX}"].to_numpy()

    diff = pred - true
    rmse = float(np.sqrt(np.mean(diff**2)))
    ss_res = float(np.sum(diff**2))
    ss_tot = float(np.sum((true - np.mean(true)) ** 2))
    r2 = float("nan") if ss_tot == 0.0 else float(1.0 - ss_res / ss_tot)

    auc_matlab = float(np.trapezoid(true, merged["time_days"]))
    auc_python = float(np.trapezoid(pred, merged["time_days"]))

    return {
        "rmse": rmse,
        "r2": r2,
        "delta_auc": float(auc_python - auc_matlab),
    }


def _tgi(frame: pd.DataFrame) -> float:
    final = frame.iloc[-1]["tumour_volume_l"]
    baseline = frame.iloc[-1]["tumour_volume_l_control"]
    if baseline <= 0:
        return float("nan")
    return 1.0 - float(final / baseline)


def _prepare_tgi_frame(matlab_ctl: pd.DataFrame, matlab_trt: pd.DataFrame) -> pd.DataFrame:
    merged = matlab_trt.merge(
        matlab_ctl[["time_days", "tumour_volume_l"]],
        on="time_days",
        suffixes=("", "_control"),
    )
    return merged


def _compute_tgi_delta(
    matlab_ctl: pd.DataFrame,
    matlab_trt: pd.DataFrame,
    python_ctl: pd.DataFrame,
    python_trt: pd.DataFrame,
) -> float:
    matlab_full = _prepare_tgi_frame(matlab_ctl, matlab_trt)
    python_full = _prepare_tgi_frame(python_ctl, python_trt)
    return (_tgi(python_full) - _tgi(matlab_full)) * 100.0


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Quantify MATLAB vs Python agreement for tutorial scenarios.")
    parser.add_argument(
        "--matlab-root",
        type=Path,
        default=Path("artifacts/matlab_validation"),
        help="Directory containing *_matlab.csv trajectories",
    )
    parser.add_argument(
        "--python-root",
        type=Path,
        default=Path("artifacts/validation"),
        help="Directory containing reference CSV outputs",
    )
    parser.add_argument(
        "--metrics-csv",
        type=Path,
        default=None,
        help="Optional destination CSV for per-observable metrics",
    )
    parser.add_argument(
        "--tgi-csv",
        type=Path,
        default=None,
        help="Optional destination CSV for delta TGI results",
    )
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)

    scenarios = [
        Scenario("example1_control", has_control=True),
        Scenario("example1_treated", has_control=False),
        Scenario("example2_treated", has_control=False),
    ]

    records: list[dict[str, object]] = []
    for scenario in scenarios:
        matlab_frame = _load_frame(args.matlab_root / f"{scenario.name}_matlab.csv")
        python_frame = _load_frame(args.python_root / f"{scenario.name}_reference.csv")

        for observable in scenario.observables:
            metrics = _compute_metrics(matlab_frame, python_frame, observable)
            records.append(
                {
                    "scenario": scenario.name,
                    "observable": observable,
                    **metrics,
                }
            )

    metrics_df = pd.DataFrame(records)
    if args.metrics_csv:
        args.metrics_csv.parent.mkdir(parents=True, exist_ok=True)
        metrics_df.to_csv(args.metrics_csv, index=False)
    else:
        print(metrics_df)

    if args.tgi_csv:
        matlab_ctl = _load_frame(args.matlab_root / "example1_control_matlab.csv")
        matlab_trt = _load_frame(args.matlab_root / "example1_treated_matlab.csv")
        python_ctl = _load_frame(args.python_root / "example1_control_reference.csv")
        python_trt = _load_frame(args.python_root / "example1_treated_reference.csv")
        delta_tgi = _compute_tgi_delta(matlab_ctl, matlab_trt, python_ctl, python_trt)
        pd.DataFrame(
            [
                {
                    "scenario": "example1_treated",
                    "delta_tgi_pct": delta_tgi,
                }
            ]
        ).to_csv(args.tgi_csv, index=False)

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
