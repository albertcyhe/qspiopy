"""Validate frozen SimBiology snapshots against MATLAB references."""

from __future__ import annotations

import argparse
import json
import logging
import math
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from scripts.compare_audits import compare_doses, compare_events
from scripts.scenario_registry import a_series, doses_to_entries, microdose
from src.offline.initial_conditions import ICOptions
from src.offline.frozen_model import EVENT_LOG_FIELDS, load_frozen_model, simulate_frozen_model


@dataclass
class Scenario:
    name: str
    parameters: Tuple[Path, ...]
    therapy: str
    snapshot: str
    stop_time: float = 400.0
    sample_interval_hours: Optional[float] = None
    custom_doses: Optional[Tuple] = None
    context_outputs: Optional[Dict[str, str]] = None


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

A1_SPEC = a_series()[0]
A1_DOSES = tuple(doses_to_entries(A1_SPEC.doses))
MICRO_SPEC = microdose()
MICRO_DOSES = tuple(doses_to_entries(MICRO_SPEC.doses))
SCENARIO_REGISTRY["A1"] = Scenario(
    "A1",
    tuple(),
    "anti_pd1",
    snapshot=A1_SPEC.snapshot,
    stop_time=A1_SPEC.days,
    sample_interval_hours=A1_SPEC.sample_interval_hours,
    custom_doses=A1_DOSES,
    context_outputs=dict(A1_SPEC.context_outputs),
)
SCENARIO_REGISTRY[MICRO_SPEC.name] = Scenario(
    MICRO_SPEC.name,
    tuple(),
    MICRO_SPEC.therapy,
    snapshot=MICRO_SPEC.snapshot,
    stop_time=MICRO_SPEC.days,
    sample_interval_hours=MICRO_SPEC.sample_interval_hours,
    custom_doses=MICRO_DOSES,
    context_outputs=dict(MICRO_SPEC.context_outputs),
)


def _parse_param_overrides(pairs: Iterable[str]) -> Dict[str, float]:
    overrides: Dict[str, float] = {}
    for pair in pairs:
        if "=" not in pair:
            raise SystemExit(f"Invalid --param-override '{pair}'; expected name=value")
        name, value = pair.split("=", 1)
        try:
            overrides[name.strip()] = float(value)
        except ValueError as exc:  # pragma: no cover - CLI validation
            raise SystemExit(f"Invalid numeric value in --param-override '{pair}': {exc}") from exc
    return overrides


def _collect_pk_invariants(model) -> Dict[str, float]:
    invariants: Dict[str, float] = {}
    params = model.parameters
    compartments = model.compartments

    def ratio(param_key: str, compartment_key: str, label: str) -> None:
        value = params.get(param_key)
        volume = compartments.get(compartment_key)
        if value is None or volume in (None, 0.0):
            return
        invariants[label] = float(value) / float(volume)

    ratio("k_cl_nivolumab", "V_C", "kel_nivolumab")
    ratios = [
        ("q_T_nivolumab", "V_C", "k12_nivolumab_T"),
        ("q_T_nivolumab", "V_T", "k21_nivolumab_T"),
        ("q_LN_nivolumab", "V_C", "k12_nivolumab_LN"),
        ("q_LN_nivolumab", "V_LN", "k21_nivolumab_LN"),
        ("q_P_nivolumab", "V_C", "k12_nivolumab_P"),
        ("q_P_nivolumab", "V_P", "k21_nivolumab_P"),
    ]
    for param_key, compartment_key, label in ratios:
        ratio(param_key, compartment_key, label)
    return invariants


def _estimate_kel(frame: pd.DataFrame, last_dose_day: Optional[float]) -> Optional[float]:
    if "drug_plasma_molar" not in frame.columns:
        return None
    times = frame["time_days"].to_numpy(dtype=float)
    conc = frame["drug_plasma_molar"].to_numpy(dtype=float)
    mask = conc > 0
    if last_dose_day is not None:
        mask &= times >= (last_dose_day + 1.0)
    tail_times = times[mask]
    tail_conc = conc[mask]
    if tail_conc.size < 4:
        return None
    log_conc = np.log(tail_conc)
    slope, _ = np.polyfit(tail_times, log_conc, 1)
    kel = -float(slope)
    if kel <= 0 or not np.isfinite(kel):
        return None
    return kel


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
    ic_mode: str,
    ic_options: ICOptions,
    module_blocks: Sequence[str],
    param_overrides: Dict[str, float],
) -> Tuple[Dict[str, Path], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_path = output_dir / f"{scenario.name}_reference.csv"
    if not reference_path.exists():
        raise SystemExit(f"Reference trajectory {reference_path} is missing; regenerate the MATLAB snapshot outputs.")
    reference_frame = pd.read_csv(reference_path, comment="#")

    events_path = output_dir / f"events_{scenario.name}_python.csv"
    surrogate_events_path = output_dir / f"{scenario.name}_events_surrogate.csv"

    model_meta = load_frozen_model(scenario.snapshot)
    if scenario.sample_interval_hours is not None:
        interval = scenario.sample_interval_hours
        time_unit = model_meta.time_unit
        if time_unit.lower().startswith("day") and math.isclose(interval, 24.0, rel_tol=0.0, abs_tol=1e-9):
            raise SystemExit(
                f"NUMERIC_GATE_FAIL scenario={scenario.name} "
                "sample_interval_hours=24 while model time unit is 'day'; use hour-based sampling."
            )

    collect_audit = collect_events
    candidate_events: Optional[List[Dict[str, object]]] = [] if (emit_diagnostics or collect_events) else None
    dose_audit: Optional[List[Dict[str, object]]] = [] if collect_audit else None
    sample_interval = scenario.sample_interval_hours
    simulate_kwargs = {}
    if sample_interval is not None:
        simulate_kwargs["sample_interval_hours"] = sample_interval
    if scenario.custom_doses is not None:
        simulate_kwargs["custom_doses"] = scenario.custom_doses
    if scenario.context_outputs is not None:
        simulate_kwargs["context_outputs"] = scenario.context_outputs
    if collect_audit:
        simulate_kwargs["dose_audit"] = dose_audit
    if collect_events:
        simulate_kwargs["event_log"] = candidate_events
    elif emit_diagnostics:
        simulate_kwargs["event_log"] = candidate_events
    simulate_kwargs["ic_mode"] = ic_mode
    if ic_mode == "target_volume":
        simulate_kwargs["ic_options"] = ic_options
    if module_blocks:
        simulate_kwargs["module_blocks"] = list(module_blocks)
    if param_overrides:
        simulate_kwargs["param_overrides"] = dict(param_overrides)
    result = simulate_frozen_model(
        scenario.snapshot,
        days=scenario.stop_time,
        therapy=scenario.therapy,
        seed=seed,
        emit_diagnostics=emit_diagnostics,
        run_label=scenario.name,
        **simulate_kwargs,
    )
    surrogate_frame = result.to_frame()
    events_df = pd.DataFrame(candidate_events or [], columns=EVENT_LOG_FIELDS)
    if not events_df.empty:
        events_df["time_hours"] = events_df["time_fire"] * 24.0
        events_df["index_in_model"] = events_df["event_index"]
        events_df.insert(0, "scenario", scenario.name)
    elif emit_diagnostics or collect_events:
        events_df = pd.DataFrame(columns=["scenario", *EVENT_LOG_FIELDS])
    if (emit_diagnostics or collect_events):
        if not events_df.empty:
            events_df.to_csv(events_path, index=False)
            events_df.to_csv(surrogate_events_path, index=False)
        else:
            if not events_path.exists():
                events_df.to_csv(events_path, index=False)
            if not surrogate_events_path.exists():
                events_df.to_csv(surrogate_events_path, index=False)

    surrogate_path = output_dir / f"{scenario.name}_surrogate.csv"
    surrogate_frame.to_csv(surrogate_path, index=False)
    if dose_audit is not None:
        audit_df = pd.DataFrame(dose_audit, columns=[
            "time_days",
            "time_hours",
            "dose_name",
            "target",
            "interpreted_dimension",
            "units",
            "compartment",
            "compartment_volume_l",
            "delta_state_value",
            "delta_amount_mol",
            "amount_moles",
            "amount_mg",
            "time_unit",
        ])
        audit_df.to_csv(output_dir / f"{scenario.name}_surrogate_dose_audit.csv", index=False)
    else:
        audit_file = output_dir / f"{scenario.name}_surrogate_dose_audit.csv"
        if not audit_file.exists():
            pd.DataFrame(columns=[
                "time_days",
                "time_hours",
                "dose_name",
                "target",
                "interpreted_dimension",
                "units",
                "compartment",
                "compartment_volume_l",
                "delta_state_value",
                "delta_amount_mol",
                "amount_moles",
                "amount_mg",
                "time_unit",
            ]).to_csv(audit_file, index=False)

    if dump_t0:
        model = model_meta
        t0_path = model.source_dir / "equations_eval_t0.csv"
        if t0_path.is_file():
            shutil.copyfile(t0_path, output_dir / f"{scenario.name}_equations_eval_t0.csv")
        else:
            logging.warning("equations_eval_t0.csv not found under %s", model.source_dir)

    grid_hours = surrogate_frame["time_days"].to_numpy(dtype=float) * 24.0
    last_dose_day: Optional[float] = None
    if scenario.custom_doses:
        last_dose_day = max(float(entry.start_time) for entry in scenario.custom_doses)

    grid_meta = {
        "scenario": scenario.name,
        "grid_first_hours": np.round(grid_hours[:6], 6).tolist(),
        "grid_unique_step_hours": np.unique(np.round(np.diff(grid_hours), 6))[:3].tolist()
        if grid_hours.size > 1
        else [],
    }
    pk_meta = {"scenario": scenario.name, **_collect_pk_invariants(model_meta)}
    pk_meta["kel_param"] = model_meta.parameters.get("k_cl_nivolumab")
    kel_sur = _estimate_kel(surrogate_frame, last_dose_day)
    kel_ref = _estimate_kel(reference_frame, last_dose_day)
    if kel_sur is not None:
        pk_meta["kel_estimate_surrogate"] = kel_sur
    if kel_ref is not None:
        pk_meta["kel_estimate_reference"] = kel_ref

    return (
        {"surrogate": surrogate_path, "reference": reference_path},
        surrogate_frame,
        reference_frame,
        events_df,
        grid_meta,
        pk_meta,
    )


def _performance_benchmark(parameters: Tuple[Path, ...], replicates: int) -> Dict[str, float]:
    start = time.perf_counter()
    simulate_frozen_model("example1", days=400.0, therapy="anti_pd1", sample_interval_hours=12.0)
    ref_time = time.perf_counter() - start

    start = time.perf_counter()
    for _ in range(replicates):
        simulate_frozen_model("example1", days=400.0, therapy="anti_pd1", sample_interval_hours=12.0)
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


def _enforce_numeric_gates(
    output_dir: Path,
    run_data: List[
        Tuple[str, pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, object], Dict[str, float]]
    ],
) -> None:
    key_columns = ["tumour_volume_l", "pd1_occupancy", "tcell_density_per_ul"]
    rel_l2_threshold = 1e-3
    max_re_threshold = 5e-3
    event_threshold_hours = 0.5

    for scenario_name, surrogate, reference, events_df, _, _ in run_data:
        issues: List[str] = []
        for column in key_columns:
            if column not in surrogate.columns or column not in reference.columns:
                issues.append(f"missing column '{column}' in outputs")
                continue
            rel_l2 = _compute_rel_l2(surrogate[column], reference[column])
            if rel_l2 > rel_l2_threshold:
                issues.append(f"column={column} rel_L2={rel_l2:.6g}>{rel_l2_threshold}")
            max_rel = _compute_max_rel(surrogate[column], reference[column])
            if max_rel > max_re_threshold:
                issues.append(f"column={column} maxRE={max_rel:.6g}>{max_re_threshold}")

        reference_events_path = output_dir / f"events_{scenario_name}_reference.csv"
        if not reference_events_path.is_file():
            alt_ref = output_dir / f"{scenario_name}_events_reference.csv"
            if alt_ref.is_file():
                reference_events_path = alt_ref
        surrogate_events_path = output_dir / f"{scenario_name}_events_surrogate.csv"
        if not surrogate_events_path.is_file():
            alt_sur = output_dir / f"events_{scenario_name}_surrogate.csv"
            if alt_sur.is_file():
                surrogate_events_path = alt_sur
        if reference_events_path.is_file() and surrogate_events_path.is_file():
            event_issues = compare_events(surrogate_events_path, reference_events_path, tol_h=event_threshold_hours)
            issues.extend(event_issues)
        else:
            issues.append("missing event audit files for comparison")

        surrogate_dose_path = output_dir / f"{scenario_name}_surrogate_dose_audit.csv"
        reference_dose_path = output_dir / f"{scenario_name}_reference_dose_audit.csv"
        if surrogate_dose_path.is_file() and reference_dose_path.is_file():
            dose_issues = compare_doses(surrogate_dose_path, reference_dose_path)
            issues.extend(dose_issues)
        else:
            issues.append("missing dose audit files for comparison")

        if issues:
            print(f"NUMERIC_GATE_FAIL scenario={scenario_name}")
            for item in issues:
                print(f" - {item}")
            raise SystemExit(3)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate the Python surrogate against the reference equations")
    parser.add_argument("--output", type=Path, default=Path("artifacts/validation"))
    parser.add_argument(
        "--benchmark-replicates",
        type=int,
        default=0,
        help="Number of runs for the optional performance benchmark (0 disables it)",
    )
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
    parser.add_argument("--ic-mode", choices=["snapshot", "target_volume"], default="target_volume")
    parser.add_argument("--ic-target-diam-cm", type=float, default=0.5)
    parser.add_argument(
        "--ic-reset-policy",
        choices=["all_zero", "cancer_only", "custom"],
        default="cancer_only",
        help="State reset policy for target_volume ICs",
    )
    parser.add_argument(
        "--ic-preserve-pattern",
        action="append",
        default=[],
        help="fnmatch pattern(s) of state identifiers to preserve when reset_policy=custom",
    )
    parser.add_argument("--ic-max-days", type=float, default=150.0)
    parser.add_argument("--ic-max-wall-seconds", type=float, default=20.0)
    parser.add_argument(
        "--module-block",
        action="append",
        default=[],
        help="Repeated-assignment targets to disable before simulation",
    )
    parser.add_argument(
        "--param-override",
        action="append",
        default=[],
        help="Parameter override in the form name=value (may be repeated)",
    )
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

    preserve_patterns = tuple(
        pattern.strip().lower()
        for pattern in (args.ic_preserve_pattern or [])
        if pattern and pattern.strip()
    )
    ic_options = ICOptions(
        target_diameter_cm=args.ic_target_diam_cm,
        reset_policy=args.ic_reset_policy,
        preserve_patterns=preserve_patterns,
        max_days=args.ic_max_days,
        max_wall_seconds=args.ic_max_wall_seconds,
    )
    param_overrides = _parse_param_overrides(args.param_override or [])
    module_blocks = [block.strip() for block in (args.module_block or []) if block and block.strip()]

    results: List[Dict[str, Path]] = []
    metrics_records: List[Dict[str, object]] = []
    worst_alignment: Optional[Dict[str, object]] = None
    observables = ["tumour_volume_l", "pd1_occupancy", "tcell_density_per_ul"]
    run_data: List[
        Tuple[str, pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, object], Dict[str, float]]
    ] = []

    for scenario in scenarios:
        paths, surrogate_frame, reference_frame, events_df, grid_meta, pk_meta = _run_scenario(
            scenario,
            args.output,
            emit_diagnostics=args.emit_diagnostics,
            seed=args.seed,
            collect_events=args.numeric_gates,
            dump_t0=args.dump_t0,
            ic_mode=args.ic_mode,
            ic_options=ic_options,
            module_blocks=module_blocks,
            param_overrides=param_overrides,
        )
        results.append({"scenario": scenario.name, **paths})
        run_data.append((scenario.name, surrogate_frame, reference_frame, events_df, grid_meta, pk_meta))

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

    example1_paths = [
        args.output / "example1_treated_surrogate.csv",
        args.output / "example1_treated_reference.csv",
        args.output / "example1_control_surrogate.csv",
        args.output / "example1_control_reference.csv",
    ]
    if all(path.exists() for path in example1_paths):
        example1_sur = pd.read_csv(example1_paths[0])
        example1_ref = pd.read_csv(example1_paths[1], comment="#")
        example1_sur_ctrl = pd.read_csv(example1_paths[2])
        example1_ref_ctrl = pd.read_csv(example1_paths[3], comment="#")

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

    grid_records = [meta for _, _, _, _, meta, _ in run_data]
    pk_records = [pk for _, _, _, _, _, pk in run_data]
    if grid_records:
        grid_frame = pd.DataFrame(grid_records)
        grid_frame.to_csv(args.output / "grid_info.csv", index=False)
    if pk_records:
        pk_frame = pd.DataFrame(pk_records)
        pk_frame.to_csv(args.output / "pk_invariants.csv", index=False)

    if args.benchmark_replicates > 0:
        benchmark = _performance_benchmark(
            (Path("parameters/example1_parameters.json"),),
            args.benchmark_replicates,
        )
        benchmark_payload = {"replicates": args.benchmark_replicates, **benchmark}
        if grid_records:
            benchmark_payload["grid_debug"] = grid_records[0]
        if pk_records:
            benchmark_payload["pk_invariants"] = pk_records[0]
        with (args.output / "performance.json").open("w", encoding="utf8") as handle:
            json.dump(benchmark_payload, handle, indent=2)

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
