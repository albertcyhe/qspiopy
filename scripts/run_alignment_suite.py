"""Run extended alignment scenarios between Python surrogate and MATLAB."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

from scripts.validate_surrogate import _compute_metrics
from src.offline.entities import BASE_HEADER
from src.offline.frozen_model import DoseEntry, simulate_frozen_model


DRUG_INFO: Dict[str, Dict[str, float | str]] = {
    "nivolumab": {"mw_mg_per_mol": 1.436e8, "target": "V_C.nivolumab"},
    "durvalumab": {"mw_mg_per_mol": 1.436e8, "target": "V_C.durvalumab"},
    "ipilimumab": {"mw_mg_per_mol": 1.486349e8, "target": "V_C.ipililumab"},
}


@dataclass(frozen=True)
class DoseSpec:
    drug: str
    start_day: float
    interval_days: float
    repeat_count: int
    amount_mg: Optional[float] = None
    amount_mg_per_kg: Optional[float] = None
    label: Optional[str] = None

    def amount_mg_weighted(self, patient_weight_kg: float) -> float:
        if self.amount_mg is not None:
            return float(self.amount_mg)
        if self.amount_mg_per_kg is None:
            raise ValueError(f"DoseSpec for {self.drug} missing amount specification")
        return float(self.amount_mg_per_kg) * float(patient_weight_kg)


@dataclass(frozen=True)
class ScenarioSpec:
    scenario_id: str
    snapshot: str
    model_script: str
    therapy: str
    stop_time_days: float
    patient_weight_kg: float
    doses: Sequence[DoseSpec]
    sample_interval_hours: float = 0.5
    extra_outputs: Dict[str, str] = field(default_factory=dict)


def _resolve_repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _build_dose_entries(spec: ScenarioSpec, *, index_offset: int = 10_000) -> List[DoseEntry]:
    entries: List[DoseEntry] = []
    for idx, dose in enumerate(spec.doses, start=1):
        drug = dose.drug.lower()
        if drug not in DRUG_INFO:
            raise ValueError(f"Unsupported drug '{dose.drug}'")
        info = DRUG_INFO[drug]
        mw = float(info["mw_mg_per_mol"])
        target = str(info["target"])
        amount_mg = dose.amount_mg_weighted(spec.patient_weight_kg)
        dose_moles = amount_mg / mw
        entry = DoseEntry(
            index=index_offset + idx,
            name=dose.label or f"{spec.scenario_id}_{drug}_{idx}",
            dose_type="RepeatDose",
            target=target,
            amount=dose_moles,
            amount_units="mole",
            start_time=float(dose.start_day),
            interval=float(dose.interval_days),
            repeat_count=int(dose.repeat_count),
            rate=None,
            rate_units="",
            duration=None,
        )
        entries.append(entry)
    return entries


def _serialize_doses_for_matlab(spec: ScenarioSpec) -> List[Dict[str, float | str]]:
    payload: List[Dict[str, float | str]] = []
    for idx, dose in enumerate(spec.doses, start=1):
        amount_mg = dose.amount_mg_weighted(spec.patient_weight_kg)
        payload.append(
            {
                "name": dose.label or f"{spec.scenario_id}_{dose.drug}_{idx}",
                "drug": dose.drug,
                "amount_mg": amount_mg,
                "start_time": float(dose.start_day),
                "interval": float(dose.interval_days),
                "repeat_count": int(dose.repeat_count),
            }
        )
    return payload


def _simulate_python(spec: ScenarioSpec, *, doses: List[DoseEntry]):
    result = simulate_frozen_model(
        spec.snapshot,
        days=spec.stop_time_days,
        therapy=spec.therapy,
        sample_interval_hours=spec.sample_interval_hours,
        custom_doses=doses,
        context_outputs=spec.extra_outputs,
    )
    return result


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


def _compute_auc(series: Iterable[float], time: Iterable[float]) -> float:
    return float(np.trapezoid(series, time))


def _scenario_registry() -> Dict[str, ScenarioSpec]:
    return {
        "A1": ScenarioSpec(
            scenario_id="A1",
            snapshot="example1",
            model_script="example1",
            therapy="anti_pd1",
            stop_time_days=84.0,
            patient_weight_kg=70.0,
            doses=[
                DoseSpec(drug="nivolumab", amount_mg=200.0, start_day=0.0, interval_days=21.0, repeat_count=3),
            ],
            extra_outputs={
                "drug_plasma_molar": "V_C.nivolumab",
                "drug_tumor_molar": "V_T.nivolumab",
                "cd8_tumor_cells": "V_T.T1",
                "treg_tumor_cells": "V_T.T0",
            },
        ),
        "A2": ScenarioSpec(
            scenario_id="A2",
            snapshot="example1",
            model_script="example1",
            therapy="anti_pd1",
            stop_time_days=84.0,
            patient_weight_kg=70.0,
            doses=[
                DoseSpec(drug="nivolumab", amount_mg=240.0, start_day=0.0, interval_days=14.0, repeat_count=5),
            ],
            extra_outputs={
                "drug_plasma_molar": "V_C.nivolumab",
                "drug_tumor_molar": "V_T.nivolumab",
                "cd8_tumor_cells": "V_T.T1",
                "treg_tumor_cells": "V_T.T0",
            },
        ),
        "A3": ScenarioSpec(
            scenario_id="A3",
            snapshot="example1",
            model_script="example1",
            therapy="anti_pd1",
            stop_time_days=84.0,
            patient_weight_kg=70.0,
            doses=[
                DoseSpec(drug="nivolumab", amount_mg=400.0, start_day=0.0, interval_days=42.0, repeat_count=1),
            ],
            extra_outputs={
                "drug_plasma_molar": "V_C.nivolumab",
                "drug_tumor_molar": "V_T.nivolumab",
                "cd8_tumor_cells": "V_T.T1",
                "treg_tumor_cells": "V_T.T0",
            },
        ),
        "A4": ScenarioSpec(
            scenario_id="A4",
            snapshot="example1",
            model_script="example1",
            therapy="anti_pd1",
            stop_time_days=84.0,
            patient_weight_kg=70.0,
            doses=[
                DoseSpec(drug="nivolumab", amount_mg_per_kg=3.0, start_day=0.0, interval_days=14.0, repeat_count=5),
            ],
            extra_outputs={
                "drug_plasma_molar": "V_C.nivolumab",
                "drug_tumor_molar": "V_T.nivolumab",
                "cd8_tumor_cells": "V_T.T1",
                "treg_tumor_cells": "V_T.T0",
            },
        ),
        "A5": ScenarioSpec(
            scenario_id="A5",
            snapshot="example1",
            model_script="example1",
            therapy="anti_pd1",
            stop_time_days=84.0,
            patient_weight_kg=70.0,
            doses=[
                DoseSpec(drug="nivolumab", amount_mg=800.0, start_day=0.0, interval_days=21.0, repeat_count=0, label="load"),
                DoseSpec(drug="nivolumab", amount_mg=200.0, start_day=21.0, interval_days=21.0, repeat_count=2, label="maintenance"),
            ],
            extra_outputs={
                "drug_plasma_molar": "V_C.nivolumab",
                "drug_tumor_molar": "V_T.nivolumab",
                "cd8_tumor_cells": "V_T.T1",
                "treg_tumor_cells": "V_T.T0",
            },
        ),
        "A6": ScenarioSpec(
            scenario_id="A6",
            snapshot="example1",
            model_script="example1",
            therapy="anti_pd1",
            stop_time_days=84.0,
            patient_weight_kg=70.0,
            doses=[
                DoseSpec(drug="nivolumab", amount_mg=50.0, start_day=0.0, interval_days=42.0, repeat_count=1),
            ],
            extra_outputs={
                "drug_plasma_molar": "V_C.nivolumab",
                "drug_tumor_molar": "V_T.nivolumab",
                "cd8_tumor_cells": "V_T.T1",
                "treg_tumor_cells": "V_T.T0",
            },
        ),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Extended Python vs MATLAB alignment runner")
    parser.add_argument("--output", type=Path, default=Path("artifacts/extended_validation"))
    parser.add_argument("--matlab-cli", type=Path, default=Path("/Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab"))
    parser.add_argument("--scenarios", nargs="+", choices=list(_scenario_registry().keys()) + ["all"])
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
        dose_entries = _build_dose_entries(spec)
        result = _simulate_python(spec, doses=dose_entries)
        surrogate_path = output_dir / f"{spec.scenario_id}_surrogate.csv"
        result.save_csv(surrogate_path, order="contract", include_header_manifest=True)
        surrogate_df = result.to_frame()

        reference_path = output_dir / f"{spec.scenario_id}_reference.csv"
        if not args.skip_matlab:
            config_payload = {
                "scenario_id": spec.scenario_id,
                "model_script": spec.model_script,
                "stop_time_days": spec.stop_time_days,
                "sample_interval_hours": spec.sample_interval_hours,
                "patient_weight_kg": spec.patient_weight_kg,
                "rtol": 1e-6,
                "atol": 1e-12,
                "max_step": None,
                "doses": _serialize_doses_for_matlab(spec),
            }
            config_path = output_dir / f"{spec.scenario_id}_scenario.json"
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
            per_column.insert(0, "scenario", spec.scenario_id)
            metrics_records.extend(per_column.to_dict(orient="records"))

            # Add aggregate metrics replicating validate_surrogate for core observables
            for observable in ["tumour_volume_l", "pd1_occupancy", "tcell_density_per_ul"]:
                met = _compute_metrics(surrogate, reference, observable)
                metrics_records.append(
                    {
                        "scenario": spec.scenario_id,
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
