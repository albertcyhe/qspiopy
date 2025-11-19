"""Multi-scenario PD-1 driver comparison utility (Step 1 of the new plan)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _latest_flat_debug(path_dir: Path, scenario: str) -> Path:
    matches = sorted(path_dir.glob(f"{scenario}_flat_debug_*.csv"))
    if not matches:
        raise FileNotFoundError(f"No flat-debug CSV found for {scenario} under {path_dir}")
    return matches[-1]


def _load_reference(path_dir: Path, scenario: str) -> pd.DataFrame:
    ref_path = path_dir / f"{scenario}_reference.csv"
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference trajectory missing: {ref_path}")
    return pd.read_csv(ref_path, comment="#")


def _load_surrogate(path_dir: Path, scenario: str) -> pd.DataFrame:
    sur_path = path_dir / f"{scenario}_surrogate.csv"
    if not sur_path.exists():
        raise FileNotFoundError(f"Surrogate trajectory missing: {sur_path}")
    return pd.read_csv(sur_path)


def _occupancy_series(frame: pd.DataFrame) -> pd.Series:
    if "pd1_occupancy" in frame.columns:
        return frame["pd1_occupancy"]
    if "pd1_occupancy_ctx" in frame.columns:
        return frame["pd1_occupancy_ctx"]
    if "H_PD1_C1" in frame.columns:
        return frame["H_PD1_C1"]
    raise KeyError("No PD-1 occupancy column found in flat debug frame")


def _extract_metrics(frame: pd.DataFrame) -> Dict[str, float]:
    pk_state = frame["pd1_alignment_pk_state"].to_numpy()
    occ = _occupancy_series(frame).to_numpy()
    pk_max = float(np.nanmax(pk_state))
    occ_max = float(np.nanmax(occ))
    occ_min = float(np.nanmin(occ))
    occ_span = occ_max - occ_min
    saturation_time = float("nan")
    mask = occ >= 0.95 * max(occ_max, 1e-12)
    if mask.any():
        saturation_time = float(frame.loc[mask, "time_days"].iloc[0])
    return {
        "pk_state_max": pk_max,
        "pk_state_span": float(np.nanmax(pk_state) - np.nanmin(pk_state)),
        "occ_max": occ_max,
        "occ_min": occ_min,
        "occ_span": occ_span,
        "first_time_95pct": saturation_time,
    }


def _plot_scenario(
    scenario: str,
    flat_frame: pd.DataFrame,
    reference: pd.DataFrame,
    surrogate: pd.DataFrame,
    *,
    output_dir: Path,
) -> Path:
    time_debug = flat_frame["time_days"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    axes[0].plot(time_debug, flat_frame["pd1_alignment_pk_state"], label="PK state", color="#1f77b4")
    axes[0].plot(
        time_debug,
        flat_frame["pd1_alignment_concentration_M"],
        label="Concentration (M)",
        color="#ff7f0e",
    )
    axes[0].set_ylabel("Driver PK")
    axes[0].grid(True, alpha=0.3, linestyle="--")
    axes[0].legend(loc="best")

    axes[1].plot(
        time_debug,
        flat_frame.get("pd1_filter_input", np.nan),
        label="Filter input",
        color="#2ca02c",
    )
    axes[1].plot(
        time_debug,
        flat_frame.get("pd1_filter_output", np.nan),
        label="Filter output",
        color="#d62728",
    )
    axes[1].plot(
        time_debug,
        flat_frame.get("pd1_filter_surface_density", np.nan),
        label="Surface density",
        color="#9467bd",
        linestyle="--",
    )
    axes[1].set_ylabel("Effective binding signal")
    axes[1].grid(True, alpha=0.3, linestyle="--")
    axes[1].legend(loc="best")

    axes[2].plot(
        surrogate["time_days"],
        surrogate["pd1_occupancy"],
        label="Surrogate",
        color="#d62728",
    )
    axes[2].plot(
        reference["time_days"],
        reference["pd1_occupancy"],
        label="MATLAB reference",
        color="#1f77b4",
        linestyle="--",
    )
    axes[2].set_ylabel("PD-1 occupancy")
    axes[2].set_xlabel("Time (days)")
    axes[2].grid(True, alpha=0.3, linestyle="--")
    axes[2].legend(loc="best")

    fig.suptitle(f"{scenario} – PD-1 driver breakdown")
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{scenario}_pd1_driver_compare.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def _summarise(flat_frame: pd.DataFrame, reference: pd.DataFrame, surrogate: pd.DataFrame) -> Dict[str, float]:
    stats = _extract_metrics(flat_frame)
    ref_peak = float(np.nanmax(reference["pd1_occupancy"]))
    ref_time_span = float(reference["time_days"].iloc[-1] - reference["time_days"].iloc[0])
    stats.update(
        {
            "reference_peak": ref_peak,
            "reference_duration": ref_time_span,
        }
    )
    return stats


def _scenario_payload(
    scenario: str,
    *,
    flat_frame: pd.DataFrame,
    reference: pd.DataFrame,
    surrogate: pd.DataFrame,
    figure_path: Path,
) -> Dict[str, object]:
    stats = _summarise(flat_frame, reference, surrogate)
    stats["figure"] = str(figure_path)
    stats["flat_debug_source"] = str(flat_frame.attrs.get("source_path", ""))
    return stats


def _annotate_source(frame: pd.DataFrame, path: Path) -> None:
    frame.attrs["source_path"] = path.resolve()


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare PD-1 driver signals across multiple scenarios")
    parser.add_argument(
        "--scenarios",
        nargs="+",
        default=["A1", "A2", "A3", "A4", "A5", "A6"],
        help="Scenario IDs to inspect",
    )
    parser.add_argument("--flat-debug-dir", type=Path, default=Path("artifacts/validation"))
    parser.add_argument("--validation-dir", type=Path, default=Path("artifacts/validation"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/dev"))
    parser.add_argument("--summary-json", type=Path, help="Optional JSON path for aggregated stats")
    args = parser.parse_args(argv)

    summary: Dict[str, object] = {}
    for scenario in args.scenarios:
        flat_path = _latest_flat_debug(args.flat_debug_dir, scenario)
        flat_frame = pd.read_csv(flat_path)
        _annotate_source(flat_frame, flat_path)
        reference = _load_reference(args.validation_dir, scenario)
        surrogate = _load_surrogate(args.validation_dir, scenario)
        fig_path = _plot_scenario(
            scenario,
            flat_frame=flat_frame,
            reference=reference,
            surrogate=surrogate,
            output_dir=args.output_dir,
        )
        summary[scenario] = _scenario_payload(
            scenario,
            flat_frame=flat_frame,
            reference=reference,
            surrogate=surrogate,
            figure_path=fig_path,
        )
        print(
            f"[pd1_compare] {scenario} pk_max={summary[scenario]['pk_state_max']:.3g} "
            f"occ_max={summary[scenario]['occ_max']:.3g} first_time_95pct={summary[scenario]['first_time_95pct']:.2f}"
        )

    if args.summary_json:
        args.summary_json.parent.mkdir(parents=True, exist_ok=True)
        args.summary_json.write_text(json.dumps(summary, indent=2), encoding="utf8")
        print(f"[pd1_compare] wrote {args.summary_json}")

    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
