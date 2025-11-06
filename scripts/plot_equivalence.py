"""Generate core plots for the semantic equivalence study."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .validate_surrogate import SCENARIO_REGISTRY
from src.offline.frozen_model import simulate_frozen_model

DEFAULT_TOLERANCES = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6)


def _ensure_non_empty(path: Path, frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        raise SystemExit(f"{path} is empty – run the validation and summary steps first.")
    return frame


def _compute_max_rel_err(reference: pd.DataFrame, candidate: pd.DataFrame, columns: Sequence[str]) -> float:
    merged = candidate.merge(reference, on="time_days", suffixes=("_py", "_ref"))
    eps = 1e-12
    worst = 0.0
    for column in columns:
        py = merged[f"{column}_py"].to_numpy()
        ref = merged[f"{column}_ref"].to_numpy()
        denom = np.maximum(np.abs(ref), eps)
        rel_err = np.abs(py - ref) / denom
        if rel_err.max(initial=0.0) > worst:
            worst = float(rel_err.max())
    return worst


def plot_max_rel_err(metrics: pd.DataFrame, output: Path) -> None:
    grouped = metrics.groupby("scenario")["max_fractional_error"].apply(list)
    data = [np.array(values) for values in grouped.values]
    if not data:
        raise SystemExit("No alignment metrics available to plot.")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.boxplot(data, labels=grouped.index, showfliers=False)
    ax.set_ylabel("Max relative error")
    ax.set_yscale("log")
    ax.set_title("Trajectory alignment (max rel. error)")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_event_residuals(event_logs: pd.DataFrame, output: Path) -> None:
    if event_logs.empty:
        raise SystemExit("Event log is empty – add at least one event-enabled snapshot.")
    fig, ax = plt.subplots(figsize=(6, 4))
    event_logs.boxplot(column="time_residual", by="scenario", ax=ax)
    ax.set_ylabel("Time residual (fire - trigger)")
    ax.set_title("Event timing residuals")
    fig.suptitle("")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_overlays(validation_dir: Path, output_dir: Path, scenarios: Iterable[str], observables: int = 3) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    missing: list[str] = []
    for scenario in scenarios:
        ref_path = validation_dir / f"{scenario}_reference.csv"
        sur_path = validation_dir / f"{scenario}_surrogate.csv"
        if not ref_path.exists() or not sur_path.exists():
            missing.append(scenario)
            continue
        ref = pd.read_csv(ref_path, comment="#")
        sur = pd.read_csv(sur_path)
        common_cols = [col for col in ref.columns if col != "time_days" and col in sur.columns]
        if not common_cols:
            continue
        for column in common_cols[:observables]:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.plot(ref["time_days"], ref[column], label="reference", linewidth=2)
            ax.plot(sur["time_days"], sur[column], label="python", linestyle="--")
            ax.set_xlabel("time (days)")
            ax.set_ylabel(column)
            ax.set_title(f"{scenario} – {column}")
            ax.legend()
            fig.tight_layout()
            fig.savefig(output_dir / f"overlay_{scenario}_{column}.png")
            plt.close(fig)
    if missing:
        raise SystemExit(f"Missing reference/surrogate CSVs for scenarios: {', '.join(missing)}")


def plot_convergence(
    validation_dir: Path,
    output_dir: Path,
    scenario_name: str,
    tolerances: Iterable[float],
    seed: int | None,
) -> None:
    scenario = SCENARIO_REGISTRY[scenario_name]
    ref_path = validation_dir / f"{scenario.name}_reference.csv"
    if not ref_path.exists():
        raise SystemExit(f"Reference trajectory {ref_path} not found.")
    reference = pd.read_csv(ref_path, comment="#")

    columns = [col for col in reference.columns if col != "time_days"]
    interval = scenario.sample_interval_hours if scenario.sample_interval_hours is not None else 12.0
    baseline = simulate_frozen_model(
        scenario.snapshot,
        days=scenario.stop_time,
        therapy=scenario.therapy,
        seed=seed,
        rtol_override=1e-8,
        atol_override=1e-12,
        sample_interval_hours=interval,
    ).to_frame()

    errors = []
    tolerances = list(tolerances)
    for tol in tolerances:
        result = simulate_frozen_model(
            scenario.snapshot,
            days=scenario.stop_time,
            therapy=scenario.therapy,
            seed=seed,
            rtol_override=tol,
            atol_override=tol * 1e-2,
            sample_interval_hours=interval,
        ).to_frame()
        err = _compute_max_rel_err(baseline, result, columns)
        errors.append(max(err, 1e-16))

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(tolerances, errors, marker="o")
    ax.set_xlabel("Relative tolerance override")
    ax.set_ylabel("max rel error vs baseline")
    ax.set_title(f"Convergence sweep – {scenario.name}")
    ax.grid(True, which="both", ls=":")
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"convergence_{scenario.name}.png")
    plt.close(fig)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate equivalence plots")
    parser.add_argument("--analysis-dir", type=Path, default=Path("artifacts/analysis"))
    parser.add_argument("--validation-dir", type=Path, default=Path("artifacts/validation"))
    parser.add_argument("--plots-dir", type=Path, default=Path("plots"))
    parser.add_argument("--scenarios", nargs="+", default=list(SCENARIO_REGISTRY.keys()))
    parser.add_argument("--convergence-scenario", default="event_suite")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--tolerances",
        type=float,
        nargs="*",
        default=DEFAULT_TOLERANCES,
        help="Relative tolerances for convergence sweep",
    )
    args = parser.parse_args(argv)

    metrics = _ensure_non_empty(
        args.analysis_dir / "alignment_metrics.csv", pd.read_csv(args.analysis_dir / "alignment_metrics.csv")
    )
    event_logs = pd.read_csv(args.analysis_dir / "event_logs.csv")

    plot_max_rel_err(metrics, args.plots_dir / "max_rel_err_boxplot.png")
    plot_event_residuals(event_logs, args.plots_dir / "event_time_residuals.png")
    plot_overlays(args.validation_dir, args.plots_dir, args.scenarios)
    plot_convergence(args.validation_dir, args.plots_dir, args.convergence_scenario, args.tolerances, args.seed)

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry
    raise SystemExit(main())
