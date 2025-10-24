"""Aggregate alignment diagnostics into publication-ready tables."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable

import pandas as pd

from .validate_surrogate import SCENARIO_REGISTRY


def _ensure_exists(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"expected file {path} is missing")


def _load_metrics(path: Path) -> pd.DataFrame:
    _ensure_exists(path)
    frame = pd.read_csv(path)
    required = {"scenario", "observable", "max_fractional_error", "relative_rmse"}
    missing = required.difference(frame.columns)
    if missing:
        raise SystemExit(f"metrics.csv missing columns: {sorted(missing)}")
    return frame


def _compute_model_scale(snapshot_root: Path, snapshot: str) -> Dict[str, int]:
    directory = snapshot_root / snapshot
    if not directory.is_dir():
        raise SystemExit(f"snapshot directory {directory} not found")

    def count_csv(name: str) -> int:
        path = directory / name
        if not path.exists():
            return 0
        frame = pd.read_csv(path)
        return int(len(frame))

    return {
        "snapshot": snapshot,
        "species": count_csv("species.csv"),
        "reactions": count_csv("reactions.csv"),
        "events": count_csv("events.csv"),
        "rules": count_csv("rules.csv"),
        "doses": count_csv("doses.csv"),
    }


def summarize_alignment_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "scenario",
        "observable",
        "max_fractional_error",
        "relative_rmse",
        "rmse",
        "r2",
        "delta_auc",
    ]
    existing = [col for col in columns if col in metrics.columns]
    return metrics.loc[:, existing].copy()


def summarize_alignment_by_scenario(metrics: pd.DataFrame) -> pd.DataFrame:
    grouped = metrics.groupby("scenario")
    summary = grouped["max_fractional_error"].max().to_frame("max_rel_err")
    summary["median_rel_rmse"] = grouped["relative_rmse"].median()
    summary["p95_rel_rmse"] = grouped["relative_rmse"].quantile(0.95)
    return summary.reset_index()


def _collect_event_logs(validation_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in sorted(validation_dir.glob("events_*_python.csv")):
        frame = pd.read_csv(path)
        if frame.empty:
            frame = pd.DataFrame(columns=["scenario", "event_index", "time_fire", "time_trigger", "delay", "type", "assignments"])
        frames.append(frame)
    if not frames:
        return pd.DataFrame(
            columns=[
                "scenario",
                "event_index",
                "time_fire",
                "time_trigger",
                "delay",
                "type",
                "assignments",
                "time_residual",
                "delay_residual",
            ]
        )
    non_empty = [frame for frame in frames if not frame.empty]
    if non_empty:
        data = pd.concat(non_empty, ignore_index=True)
    else:
        columns = [
            "scenario",
            "event_index",
            "time_fire",
            "time_trigger",
            "delay",
            "type",
            "assignments",
            "time_residual",
            "delay_residual",
        ]
        return pd.DataFrame(columns=columns)
    if {"time_fire", "time_trigger", "delay"}.issubset(data.columns):
        data["time_residual"] = data["time_fire"].astype(float) - data["time_trigger"].astype(float)
        data["delay_residual"] = data["time_residual"] - data["delay"].astype(float)
    else:
        data["time_residual"] = pd.NA
        data["delay_residual"] = pd.NA
    return data


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Summarize equivalence diagnostics")
    parser.add_argument(
        "--validation-dir",
        type=Path,
        default=Path("artifacts/validation"),
        help="Directory containing validate_surrogate outputs",
    )
    parser.add_argument(
        "--snapshot-root",
        type=Path,
        default=Path("artifacts/matlab_frozen_model"),
        help="Root directory containing frozen snapshots",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/analysis"),
        help="Directory to write aggregated tables",
    )
    args = parser.parse_args(argv)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    metrics = _load_metrics(args.validation_dir / "metrics.csv")
    alignment_metrics = summarize_alignment_metrics(metrics)
    alignment_metrics.to_csv(args.output_dir / "alignment_metrics.csv", index=False)

    scenario_summary = summarize_alignment_by_scenario(metrics)
    scenario_summary.to_csv(args.output_dir / "alignment_summary.csv", index=False)

    model_rows = []
    seen_snapshots = set()
    for scenario in SCENARIO_REGISTRY.values():
        if scenario.snapshot in seen_snapshots:
            continue
        seen_snapshots.add(scenario.snapshot)
        model_rows.append(_compute_model_scale(args.snapshot_root, scenario.snapshot))
    model_scale = pd.DataFrame(model_rows)
    model_scale.to_csv(args.output_dir / "model_scale.csv", index=False)

    event_logs = _collect_event_logs(args.validation_dir)
    event_logs.to_csv(args.output_dir / "event_logs.csv", index=False)
    event_summary_path = args.output_dir / "event_summary.csv"
    if not event_logs.empty:
        event_summary = (
            event_logs.groupby("scenario")["time_residual"]
            .agg(max_time_residual="max", median_time_residual="median")
            .reset_index()
        )
        delay_summary = (
            event_logs.groupby("scenario")["delay_residual"]
            .agg(max_delay_residual=lambda s: s.abs().max(), median_delay_residual=lambda s: s.abs().median())
            .reset_index()
        )
        event_summary = event_summary.merge(delay_summary, on="scenario", how="left")
        event_summary.to_csv(event_summary_path, index=False)
    else:
        empty_summary = pd.DataFrame(
            columns=[
                "scenario",
                "max_time_residual",
                "median_time_residual",
                "max_delay_residual",
                "median_delay_residual",
            ]
        )
        empty_summary.to_csv(event_summary_path, index=False)

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
