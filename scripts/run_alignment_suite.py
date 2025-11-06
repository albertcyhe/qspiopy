"""Run extended alignment scenarios between Python surrogate and MATLAB."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from scripts.scenario_registry import ScenarioSpec, a_series, b_series, doses_to_entries
from scripts.validate_surrogate import _compute_metrics
from src.offline.entities import BASE_HEADER
from src.offline.frozen_model import simulate_frozen_model


def _resolve_repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _serialize_doses_for_matlab(spec: ScenarioSpec) -> List[Dict[str, float | str]]:
    payload: List[Dict[str, float | str]] = []
    for idx, dose in enumerate(spec.doses, start=1):
        payload.append(
            {
                "name": f"{spec.name}_dose_{idx}",
                "drug": dose.drug,
                "amount_mg": dose.amount_mg,
                "start_time": dose.time_hours / 24.0,
                "interval": dose.interval_hours / 24.0,
                "repeat_count": int(dose.repeat),
            }
        )
    return payload


def _simulate_python(spec: ScenarioSpec):
    dose_entries = doses_to_entries(spec.doses)
    audit_rows: List[Dict[str, object]] = []
    result = simulate_frozen_model(
        spec.snapshot,
        days=spec.days,
        therapy=spec.therapy,
        sample_interval_hours=spec.sample_interval_hours,
        custom_doses=dose_entries,
        context_outputs=spec.context_outputs,
        dose_audit=audit_rows,
    )
    return result, audit_rows


def _call_matlab(matlab_cli: Path, repo_root: Path, config_path: Path, output_path: Path) -> None:
    cfg_str = config_path.resolve().as_posix().replace("'", "''")
    out_str = output_path.resolve().as_posix().replace("'", "''")
    matlab_root = (repo_root / "matlab").as_posix()
    command = (
        f"cd('{matlab_root}'); "
        "addpath('.'); "
        "addpath('model'); "
        "addpath('utils'); "
        "addpath('scripts'); "
        f"run_qspio_scenario('{cfg_str}', '{out_str}');"
    )
    completed = subprocess.run(
        [str(matlab_cli), "-batch", command],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "MATLAB simulation failed with exit code "
            f"{completed.returncode}\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )


def _load_frame(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    header_path = path.with_suffix(path.suffix + ".header.txt")
    if header_path.exists():
        ordered = [line.strip() for line in header_path.read_text(encoding="utf8").splitlines() if line.strip()]
        ordered = [col for col in ordered if col in frame.columns]
        remainder = [col for col in frame.columns if col not in ordered]
        frame = frame[ordered + remainder]
    return frame


def _compare_frames(
    surrogate: pd.DataFrame,
    reference: pd.DataFrame,
    columns: Sequence[str],
    *,
    time_col: str = "time_days",
) -> pd.DataFrame:
    merged = surrogate.merge(reference, on=time_col, suffixes=("_py", "_ref"))
    records: List[Dict[str, float | str]] = []
    for column in columns:
        col_py = f"{column}_py"
        col_ref = f"{column}_ref"
        if col_py not in merged.columns or col_ref not in merged.columns:
            continue
        diff = merged[col_py].to_numpy() - merged[col_ref].to_numpy()
        ref = merged[col_ref].to_numpy()
        denom = np.linalg.norm(ref) or 1e-12
        rel_l2 = float(np.linalg.norm(diff) / denom)
        max_rel = float(np.max(np.abs(diff) / (np.abs(ref) + 1e-12)))
        records.append(
            {
                "observable": column,
                "relative_l2": rel_l2,
                "max_relative_error": max_rel,
            }
        )
    return pd.DataFrame(records)


def _scenario_registry() -> Dict[str, ScenarioSpec]:
    specs = a_series() + b_series()
    registry: Dict[str, ScenarioSpec] = {}
    for spec in specs:
        registry[spec.name] = spec
    return registry


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Extended Python vs MATLAB alignment runner")
    parser.add_argument("--output", type=Path, default=Path("artifacts/extended_validation"))
    parser.add_argument("--matlab-cli", type=Path, default=Path("/Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab"))
    available = list(_scenario_registry().keys())
    parser.add_argument("--scenarios", nargs="+", choices=available + ["all"])
    parser.add_argument("--skip-matlab", action="store_true", help="Run only Python surrogate (no MATLAB reference)")
    parser.add_argument("--emit-metrics", action="store_true", help="Write per-observable alignment metrics")
    args = parser.parse_args(argv)

    registry = _scenario_registry()
    if not args.scenarios or "all" in args.scenarios:
        scenario_ids = list(registry.keys())
    else:
        scenario_ids = args.scenarios

    output_dir = args.output
    output_dir.mkdir(parents=True, exist_ok=True)
    repo_root = _resolve_repo_root()

    metrics_records: List[Dict[str, float | str]] = []

    for scenario_id in scenario_ids:
        spec = registry[scenario_id]
        result, audit_rows = _simulate_python(spec)
        surrogate_path = output_dir / f"{spec.name}_surrogate.csv"
        result.save_csv(surrogate_path, order="contract", include_header_manifest=True)
        surrogate_df = result.to_frame()
        if audit_rows:
            audit_path = output_dir / f"{spec.name}_dose_audit.csv"
            pd.DataFrame(audit_rows).to_csv(audit_path, index=False)

        reference_path = output_dir / f"{spec.name}_reference.csv"
        if not args.skip_matlab:
            config_payload = {
                "scenario_id": spec.name,
                "label": spec.label,
                "model_script": spec.matlab_script,
                "stop_time_days": spec.days,
                "sample_interval_hours": spec.sample_interval_hours,
                "rtol": 1e-6,
                "atol": 1e-12,
                "max_step": None,
                "doses": _serialize_doses_for_matlab(spec),
            }
            config_path = output_dir / f"{spec.name}_scenario.json"
            config_path.write_text(json.dumps(config_payload, indent=2), encoding="utf8")
            _call_matlab(args.matlab_cli, repo_root, config_path, reference_path)

        if reference_path.is_file():
            surrogate = surrogate_df.copy()
            reference = _load_frame(reference_path)
            header = result.column_order()
            for column in header:
                if column not in reference.columns:
                    reference[column] = np.nan
            reference = reference[list(header) + [col for col in reference.columns if col not in header]]
            columns = [
                "cancer_cells",
                "dead_cells",
                "t_cells",
                "tumour_volume_l",
                "tumour_diameter_cm",
                "pd1_occupancy",
                "tcell_density_per_ul",
            ] + [column for column in header if column not in BASE_HEADER]
            per_column = _compare_frames(surrogate, reference, columns)
            per_column.insert(0, "scenario", spec.name)
            metrics_records.extend(per_column.to_dict(orient="records"))

            for observable in ["tumour_volume_l", "pd1_occupancy", "tcell_density_per_ul"]:
                met = _compute_metrics(surrogate, reference, observable)
                metrics_records.append(
                    {
                        "scenario": spec.name,
                        "observable": observable,
                        **met,
                    }
                )

    if args.emit_metrics and metrics_records:
        metrics_df = pd.DataFrame(metrics_records)
        metrics_df.to_csv(output_dir / "alignment_metrics_extended.csv", index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
