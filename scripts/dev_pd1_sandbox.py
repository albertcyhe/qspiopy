"""PD-1 alignment sandbox for A-series calibration (handoff step 4.3)."""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path
from typing import Iterable, List, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts.validate_surrogate import (
    SCENARIO_REGISTRY,
    Scenario,
    _compute_metrics,
    _parse_param_overrides,
    _run_scenario,
)
from src.offline.initial_conditions import ICOptions

CONTEXT_COLUMNS: Sequence[str] = (
    "pd1_alignment_pk_state",
    "pd1_alignment_concentration_M",
    "pd1_alignment_volume_l",
    "aPD1",
    "pd1_occupancy",
    "H_PD1_C1",
)


def _resolve_scenario(name: str, module_blocks: Sequence[str] | None) -> Scenario:
    try:
        scenario = SCENARIO_REGISTRY[name]
    except KeyError as exc:  # pragma: no cover - CLI validation
        raise SystemExit(f"Unknown scenario '{name}'. Available: {', '.join(sorted(SCENARIO_REGISTRY))}") from exc
    if module_blocks:
        scenario = replace(scenario, module_blocks=tuple(module_blocks))
    return scenario


def _build_context_frame(result, columns: Sequence[str]) -> pd.DataFrame:
    contexts = result.raw_contexts or ()
    if not contexts:
        return pd.DataFrame(columns=["time_days", *columns])
    rows: List[dict[str, float]] = []
    for time_val, ctx in zip(result.time_days, contexts):
        row = {"time_days": float(time_val)}
        for column in columns:
            value = ctx.get(column)
            row[column] = float(value) if value is not None else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def _plot_overlays(
    *,
    context_frame: pd.DataFrame,
    reference_frame: pd.DataFrame,
    surrogate_frame: pd.DataFrame,
    output_path: Path,
) -> None:
    plots = [
        ("pd1_occupancy", "PD-1 occupancy", True),
        ("pd1_alignment_pk_state", "PK surrogate state (a.u.)", False),
        ("pd1_alignment_concentration_M", "PK concentration surrogate (M)", False),
        ("pd1_alignment_volume_l", "Alignment tumour volume (L)", False),
    ]
    rows = len(plots)
    fig, axes = plt.subplots(rows, 1, figsize=(10, 3 * rows), sharex=True)
    if isinstance(axes, np.ndarray):
        axes = list(axes.ravel())
    elif not isinstance(axes, (list, tuple)):
        axes = [axes]

    time_ref = reference_frame["time_days"]
    time_sur = surrogate_frame["time_days"]
    time_ctx = context_frame["time_days"] if not context_frame.empty else None

    for ax, (column, title, has_reference) in zip(axes, plots):
        ax.set_title(title)
        ax.grid(True, alpha=0.3, linestyle="--", linewidth=0.6)
        if has_reference:
            if column in reference_frame.columns:
                ax.plot(time_ref, reference_frame[column], label="MATLAB reference", color="#1f77b4")
            if column in surrogate_frame.columns:
                ax.plot(time_sur, surrogate_frame[column], label="Python surrogate", color="#d62728", linestyle="--")
        else:
            if time_ctx is not None and column in context_frame.columns:
                ax.plot(time_ctx, context_frame[column], label="alignment_driver_block", color="#9467bd")
            else:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center")
        ax.legend(loc="best")

    axes[-1].set_xlabel("Time (days)")
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def _write_summary(path: Path, *, metrics: dict[str, float], outputs: dict[str, str], overrides: dict[str, float]) -> None:
    payload = {
        "metrics": metrics,
        "outputs": outputs,
        "param_overrides": overrides,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf8")


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Plot PD-1 occupancy / driver internals for a single scenario")
    parser.add_argument("--scenario", default="A1", choices=sorted(SCENARIO_REGISTRY.keys()))
    parser.add_argument("--output", type=Path, default=Path("artifacts/dev"))
    parser.add_argument(
        "--validation-dir",
        type=Path,
        default=Path("artifacts/validation"),
        help="Directory that contains <scenario>_reference.csv and receives surrogate outputs",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--module-block", action="append", default=[], help="Override module blocks for the scenario")
    parser.add_argument(
        "--param-override",
        action="append",
        default=[],
        metavar="name=value or @file.json",
        help="Optional parameter overrides passed to the simulation",
    )
    parser.add_argument("--ic-mode", choices=["snapshot", "target_volume"], default="snapshot")
    parser.add_argument("--ic-target-diam-cm", type=float, default=0.5)
    parser.add_argument("--ic-max-days", type=float, default=150.0)
    parser.add_argument("--ic-max-wall-seconds", type=float, default=20.0)

    args = parser.parse_args(argv)

    scenario = _resolve_scenario(args.scenario, args.module_block or None)
    overrides = _parse_param_overrides(args.param_override or [])
    ic_options = ICOptions(
        target_diameter_cm=args.ic_target_diam_cm,
        max_days=args.ic_max_days,
        max_wall_seconds=args.ic_max_wall_seconds,
    )

    args.output.mkdir(parents=True, exist_ok=True)
    args.validation_dir.mkdir(parents=True, exist_ok=True)

    (
        _paths,
        surrogate_frame,
        reference_frame,
        _events_df,
        _grid_meta,
        _pk_meta,
        scenario_result,
    ) = _run_scenario(
        scenario,
        args.validation_dir,
        emit_diagnostics=False,
        seed=args.seed,
        collect_events=False,
        dump_t0=False,
        ic_mode=args.ic_mode,
        ic_options=ic_options,
        module_blocks=(),
        param_overrides=overrides,
        capture_contexts=True,
    )

    metrics = _compute_metrics(surrogate_frame, reference_frame, "pd1_occupancy")

    context_frame = _build_context_frame(scenario_result, CONTEXT_COLUMNS)
    context_path = args.output / f"{scenario.name}_pd1_context.csv"
    context_frame.to_csv(context_path, index=False)

    figure_path = args.output / f"{scenario.name}_pd1_overlay.png"
    _plot_overlays(
        context_frame=context_frame,
        reference_frame=reference_frame,
        surrogate_frame=surrogate_frame,
        output_path=figure_path,
    )

    summary_path = args.output / f"{scenario.name}_pd1_summary.json"
    _write_summary(
        summary_path,
        metrics=metrics,
        outputs={
            "context_csv": str(context_path),
            "figure": str(figure_path),
            "surrogate_csv": str((args.validation_dir / f"{scenario.name}_surrogate.csv")),
            "reference_csv": str((args.validation_dir / f"{scenario.name}_reference.csv")),
        },
        overrides=overrides,
    )

    print(
        f"[pd1_sandbox] scenario={scenario.name} rel_RMSE={metrics['relative_rmse']:.3f} "
        f"max_frac_err={metrics['max_fractional_error']:.3f}"
    )
    print(f"[pd1_sandbox] figure={figure_path}")
    print(f"[pd1_sandbox] context_csv={context_path}")
    print(f"[pd1_sandbox] summary={summary_path}")
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry
    raise SystemExit(main())
