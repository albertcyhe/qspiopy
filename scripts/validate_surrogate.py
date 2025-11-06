"""Validate frozen SimBiology snapshots against MATLAB references."""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from src.offline.frozen_model import EVENT_LOG_FIELDS, load_frozen_model, simulate_frozen_model


@dataclass
class Scenario:
    name: str
    parameters: Tuple[Path, ...]
    therapy: str
    snapshot: str
    stop_time: float = 400.0


SCENARIO_REGISTRY: Dict[str, Scenario] = {
    "example1_control": Scenario(
        "example1_control",
        (Path("parameters/example1_parameters.json"),),
        "none",
        snapshot="example1",
    ),
    "example1_treated": Scenario(
        "example1_treated",
        (Path("parameters/example1_parameters.json"),),
        "anti_pd1",
        snapshot="example1",
    ),
    "example2_treated": Scenario(
        "example2_treated",
        (Path("parameters/example2_parameters.json"),),
        "anti_pd1",
        snapshot="example2",
    ),
    "event_suite": Scenario(
        "event_suite",
        (Path("parameters/event_suite_parameters.json"),),
        "none",
        snapshot="event_suite",
        stop_time=10.0,
    ),
}


def _compute_metrics(
    surrogate: pd.DataFrame,
    reference: pd.DataFrame,
    observable: str,
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


def _run_scenario(
    scenario: Scenario,
    output_dir: Path,
    *,
    emit_diagnostics: bool,
    seed: Optional[int],
    collect_events: bool,
    dump_t0: bool,
) -> Tuple[Dict[str, Path], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_path = output_dir / f"{scenario.name}_reference.csv"
    if not reference_path.exists():
        raise SystemExit(f"Reference trajectory {reference_path} is missing; regenerate the MATLAB snapshot outputs.")
    reference_frame = pd.read_csv(reference_path, comment="#")

    events_path = output_dir / f"events_{scenario.name}_python.csv"

    candidate_events: Optional[List[Dict[str, object]]] = [] if (emit_diagnostics or collect_events) else None
    result = simulate_frozen_model(
        scenario.snapshot,
        days=scenario.stop_time,
        therapy=scenario.therapy,
        seed=seed,
        emit_diagnostics=emit_diagnostics,
        run_label=scenario.name,
        event_log=candidate_events,
    )
    surrogate_frame = result.to_frame()
    events_df = pd.DataFrame(candidate_events or [], columns=EVENT_LOG_FIELDS)
    if not events_df.empty:
        events_df.insert(0, "scenario", scenario.name)
    elif emit_diagnostics or collect_events:
        events_df = pd.DataFrame(columns=["scenario", *EVENT_LOG_FIELDS])
    if (emit_diagnostics or collect_events) and not events_df.empty:
        events_df.to_csv(events_path, index=False)
    elif (emit_diagnostics or collect_events) and not events_path.exists():
        events_df.to_csv(events_path, index=False)

    surrogate_path = output_dir / f"{scenario.name}_surrogate.csv"
    surrogate_frame.to_csv(surrogate_path, index=False)

    if dump_t0:
        model = load_frozen_model(scenario.snapshot)
        t0_path = model.source_dir / "equations_eval_t0.csv"
        if t0_path.is_file():
            shutil.copyfile(t0_path, output_dir / f"{scenario.name}_equations_eval_t0.csv")
        else:
            logging.warning("equations_eval_t0.csv not found under %s", model.source_dir)

    return {"surrogate": surrogate_path, "reference": reference_path}, surrogate_frame, reference_frame, events_df


def _performance_benchmark(parameters: Tuple[Path, ...], replicates: int) -> Dict[str, float]:
    start = time.perf_counter()
    simulate_frozen_model("example1", days=400.0, therapy="anti_pd1")
    ref_time = time.perf_counter() - start

    start = time.perf_counter()
    for _ in range(replicates):
        simulate_frozen_model("example1", days=400.0, therapy="anti_pd1")
    sur_time = time.perf_counter() - start

    return {
        "reference_seconds": ref_time,
        "reference_replicates": 1,
        "surrogate_seconds": sur_time,
        "surrogate_replicates": replicates,
    }


def summarize_alignment_errors(
    frame: pd.DataFrame,
    reference: pd.DataFrame,
    columns: Iterable[str],
    *,
    time_col: str = "time_days",
    topk: int = 3,
):
    merged = frame.merge(reference, on=time_col, suffixes=("_py", "_ref"))
    rows: List[Tuple[str, float]] = []
    worst: Optional[Dict[str, float]] = None
    if merged.empty:
        return pd.DataFrame(columns=["observable", "max_rel_err"]).head(0), None
    for column in columns:
        py_values = merged[f"{column}_py"].to_numpy()
        ref_values = merged[f"{column}_ref"].to_numpy()
        denom = np.maximum(np.abs(ref_values), 1e-12)
        rel_errors = np.abs(py_values - ref_values) / denom
        idx = int(rel_errors.argmax())
        max_err = float(rel_errors[idx])
        rows.append((column, max_err))
        if worst is None or max_err > worst["max_rel_err"]:
            worst = {
                "observable": column,
                "max_rel_err": max_err,
                "time": float(merged[time_col].iloc[idx]),
                "value_py": float(py_values[idx]),
                "value_ref": float(ref_values[idx]),
            }
    table = pd.DataFrame(rows, columns=["observable", "max_rel_err"])
    table.sort_values("max_rel_err", ascending=False, inplace=True)
    return table.head(topk), worst


def _compute_rel_l2(surrogate: pd.Series, reference: pd.Series) -> float:
    diff = surrogate.to_numpy() - reference.to_numpy()
    denom = np.linalg.norm(reference.to_numpy()) or 1e-12
    return float(np.linalg.norm(diff) / denom)


def _compute_max_rel(surrogate: pd.Series, reference: pd.Series) -> float:
    diff = surrogate.to_numpy() - reference.to_numpy()
    denom = np.maximum(np.abs(reference.to_numpy()), 1e-12)
    return float(np.max(np.abs(diff) / denom))


def _enforce_numeric_gates(output_dir: Path, run_data: List[Tuple[str, pd.DataFrame, pd.DataFrame, pd.DataFrame]]) -> None:
    key_columns = ["tumour_volume_l", "pd1_occupancy", "tcell_density_per_ul"]
    rel_l2_threshold = 1e-3
    max_re_threshold = 5e-3
    event_threshold_hours = 0.5

    for scenario_name, surrogate, reference, events_df in run_data:
        for column in key_columns:
            if column not in surrogate.columns or column not in reference.columns:
                raise SystemExit(f"NUMERIC_GATE_FAIL scenario={scenario_name} missing column '{column}' in outputs")
            rel_l2 = _compute_rel_l2(surrogate[column], reference[column])
            if rel_l2 > rel_l2_threshold:
                raise SystemExit(
                    f"NUMERIC_GATE_FAIL scenario={scenario_name} column={column} rel_L2={rel_l2:.6g}>{rel_l2_threshold}"
                )
            max_rel = _compute_max_rel(surrogate[column], reference[column])
            if max_rel > max_re_threshold:
                raise SystemExit(
                    f"NUMERIC_GATE_FAIL scenario={scenario_name} column={column} maxRE={max_rel:.6g}>{max_re_threshold}"
                )

        reference_events_path = output_dir / f"events_{scenario_name}_reference.csv"
        if not reference_events_path.is_file():
            raise SystemExit(f"NUMERIC_GATE_FAIL scenario={scenario_name} missing reference event log {reference_events_path}")
        reference_events = pd.read_csv(reference_events_path)
        if len(reference_events) != len(events_df):
            raise SystemExit(
                f"NUMERIC_GATE_FAIL scenario={scenario_name} event count mismatch (reference={len(reference_events)} vs python={len(events_df)})"
            )
        if len(reference_events) > 0:
            ref_sorted = reference_events.sort_values(["event_index", "time_fire"]).reset_index(drop=True)
            cand_sorted = events_df.sort_values(["event_index", "time_fire"]).reset_index(drop=True)
            delta_hours = np.abs(cand_sorted["time_fire"] - ref_sorted["time_fire"]) * 24.0
            if (delta_hours > event_threshold_hours).any():
                worst = float(delta_hours.max())
                raise SystemExit(
                    f"NUMERIC_GATE_FAIL scenario={scenario_name} event_time_delta_hours={worst:.6g}>{event_threshold_hours}"
                )


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate the Python surrogate against the reference equations")
    parser.add_argument("--output", type=Path, default=Path("artifacts/validation"))
    parser.add_argument("--replicates", type=int, default=50, help="Number of runs for the performance benchmark")
    parser.add_argument(
        "--max-rel-err",
        type=float,
        default=1e-6,
        help="Maximum tolerated fractional error between reference and surrogate trajectories.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Deterministic seed propagated to simulations")
    parser.add_argument("--emit-diagnostics", action="store_true", help="Emit solver/banner diagnostics and error tables")
    parser.add_argument("--numeric-gates", action="store_true", help="Enforce rel_L2/maxRE/event timing thresholds")
    parser.add_argument("--dump-t0", action="store_true", help="Copy equations_eval_t0.csv into the output directory")
    scenario_choices = sorted(SCENARIO_REGISTRY.keys())
    parser.add_argument(
        "--scenarios",
        nargs="+",
        choices=scenario_choices + ["all"],
        help="Subset of scenarios to validate (default: all)",
    )
    args = parser.parse_args(argv)

    if args.emit_diagnostics:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(message)s")
    logger = logging.getLogger("validate_surrogate")

    if args.scenarios:
        if "all" in args.scenarios:
            scenario_names = scenario_choices
        else:
            scenario_names = list(dict.fromkeys(args.scenarios))
    else:
        scenario_names = scenario_choices
    scenarios = [SCENARIO_REGISTRY[name] for name in scenario_names]

    results: List[Dict[str, Path]] = []
    metrics_records: List[Dict[str, object]] = []
    worst_alignment: Optional[Dict[str, object]] = None
    observables = ["tumour_volume_l", "pd1_occupancy", "tcell_density_per_ul"]
    run_data: List[Tuple[str, pd.DataFrame, pd.DataFrame, pd.DataFrame]] = []

    for scenario in scenarios:
        paths, surrogate_frame, reference_frame, events_df = _run_scenario(
            scenario,
            args.output,
            emit_diagnostics=args.emit_diagnostics,
            seed=args.seed,
            collect_events=args.numeric_gates,
            dump_t0=args.dump_t0,
        )
        results.append({"scenario": scenario.name, **paths})
        run_data.append((scenario.name, surrogate_frame, reference_frame, events_df))

        surrogate = surrogate_frame
        reference = reference_frame

        top_table, worst_sur = summarize_alignment_errors(surrogate, reference, observables)
        if args.emit_diagnostics and not top_table.empty:
            logger.info(
                "scenario=%s surrogate_vs_reference:\n%s",
                scenario.name,
                top_table.to_string(index=False),
            )

        for observable in observables:
            metrics = _compute_metrics(surrogate, reference, observable)
            metrics_records.append(
                {
                    "scenario": scenario.name,
                    "observable": observable,
                    **metrics,
                }
            )

        if worst_sur:
            candidate = {"scenario": scenario.name, **worst_sur, "source": "candidate"}
            if worst_alignment is None or candidate["max_rel_err"] > worst_alignment["max_rel_err"]:
                worst_alignment = candidate

    example1_sur = pd.read_csv(args.output / "example1_treated_surrogate.csv")
    example1_ref = pd.read_csv(args.output / "example1_treated_reference.csv", comment="#")
    example1_sur_ctrl = pd.read_csv(args.output / "example1_control_surrogate.csv")
    example1_ref_ctrl = pd.read_csv(args.output / "example1_control_reference.csv", comment="#")

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

    if args.numeric_gates:
        _enforce_numeric_gates(args.output, run_data)

    if worst_alignment is None:
        logger.error("Validation failed: no alignment data available.")
        raise SystemExit(3)
    if worst_alignment["max_rel_err"] > args.max_rel_err:
        logger.error(
            "ALIGNMENT_FAIL source=%s scenario=%s observable=%s time=%.6g rel_err=%.3e>%.3e py=%.6g ref=%.6g",
            worst_alignment.get("source", "candidate"),
            worst_alignment["scenario"],
            worst_alignment["observable"],
            worst_alignment["time"],
            worst_alignment["max_rel_err"],
            args.max_rel_err,
            worst_alignment["value_py"],
            worst_alignment["value_ref"],
        )
        if args.emit_diagnostics:
            logger.info("Diagnostics saved under %s", args.output)
        raise SystemExit(3)

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry
    raise SystemExit(main())
