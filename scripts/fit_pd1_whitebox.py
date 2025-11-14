#!/usr/bin/env python

"""Fit PD-1 white-box parameters against MATLAB-generated training data."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import sys
import math

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.offline.modules.pd1_whitebox import PD1WhiteboxModel


@dataclass
class ScenarioData:
    name: str
    dt: np.ndarray
    concentration: np.ndarray
    target: np.ndarray
    init_context: Dict[str, float]


def load_parameters(path: Path) -> Dict[str, float]:
    payload = json.loads(path.read_text())
    params = {}
    for entry in payload:
        if not isinstance(entry, dict):
            continue
        name = entry.get("name")
        if name is None:
            continue
        value = entry.get("value")
        if value is not None:
            params[name] = float(value)
    # Fill derived PD-1 constants if missing
    d_syn = params.get("d_syn", 1.0)
    if "kon_PD1_PDL1" not in params and "k_PD1_PDL1" in params:
        params["kon_PD1_PDL1"] = params["k_PD1_PDL1"] / d_syn
    if "kon_PD1_PDL2" not in params and "k_PD1_PDL2" in params:
        params["kon_PD1_PDL2"] = params["k_PD1_PDL2"] / d_syn
    if "koff_PD1_PDL1" not in params and "k_PD1_PDL1" in params and "kd_PD1_PDL1" in params:
        params["koff_PD1_PDL1"] = params["k_PD1_PDL1"] * params["kd_PD1_PDL1"]
    if "koff_PD1_PDL2" not in params and "k_PD1_PDL2" in params and "kd_PD1_PDL2" in params:
        params["koff_PD1_PDL2"] = params["k_PD1_PDL2"] * params["kd_PD1_PDL2"]
    if "koff_PD1_aPD1" not in params and "kon_PD1_aPD1" in params and "kd_PD1_aPD1" in params:
        params["koff_PD1_aPD1"] = params["kon_PD1_aPD1"] * params["kd_PD1_aPD1"]
    if "A_Tcell" not in params and "D_Tcell" in params:
        d_t = params["D_Tcell"]
        params["A_Tcell"] = 4.0 * math.pi * (d_t / 2.0) ** 2
    if "A_cell" not in params and "D_cell" in params:
        d_c = params["D_cell"]
        params["A_cell"] = 4.0 * math.pi * (d_c / 2.0) ** 2
    return params


def prepare_training(frame: pd.DataFrame) -> List[ScenarioData]:
    scenarios: List[ScenarioData] = []
    grouped = frame.groupby("scenario_id", sort=False)
    for name, group in grouped:
        group = group.sort_values("time_days")
        times = group["time_days"].to_numpy(dtype=float)
        conc = group["drug_tumor_molar"].to_numpy(dtype=float)
        target = group["pd1_inhibition"].to_numpy(dtype=float)
        if conc.size < 2:
            continue
        dt = np.empty_like(times)
        dt[0] = 0.0
        dt[1:] = np.diff(times)
        init_context = {
            "syn_pd1_pdl1": float(group["syn_pd1_pdl1"].iloc[0]),
            "syn_pd1_pdl2": float(group["syn_pd1_pdl2"].iloc[0]),
            "syn_pd1_apd1": float(group["syn_pd1_apd1"].iloc[0]),
            "syn_pd1_apd1_pd1": float(group["syn_pd1_apd1_pd1"].iloc[0]),
        }
        scenarios.append(
            ScenarioData(
                name=str(name),
                dt=dt,
                concentration=conc,
                target=target,
                init_context=init_context,
            )
        )
    return scenarios


def simulate_scenario(
    scenario: ScenarioData,
    parameters: Dict[str, float],
) -> np.ndarray:
    context = dict(scenario.init_context)
    model = PD1WhiteboxModel.from_context(parameters, context)
    preds = np.empty_like(scenario.target)
    for idx, (dt, conc) in enumerate(zip(scenario.dt, scenario.concentration)):
        outputs = model.step(conc, max(dt, 0.0))
        preds[idx] = outputs.occupancy
    return preds


def filter_valid_scenarios(
    scenarios: Iterable[ScenarioData],
    parameters: Dict[str, float],
) -> List[ScenarioData]:
    valid: List[ScenarioData] = []
    skipped: List[str] = []
    for scenario in scenarios:
        preds = simulate_scenario(scenario, parameters)
        if np.all(np.isfinite(preds)):
            valid.append(scenario)
        else:
            skipped.append(scenario.name)
    if skipped:
        print(f"[fit_pd1] Skipping {len(skipped)} scenarios with non-finite responses: {skipped[:5]}{'...' if len(skipped) > 5 else ''}")
    return valid


def fit_parameters(
    scenarios: Iterable[ScenarioData],
    base_params: Dict[str, float],
    *,
    output_json: Path,
) -> Dict[str, float]:
    scenario_list = list(scenarios)
    total_points = sum(len(s.target) for s in scenario_list)

    def make_params(vec: np.ndarray) -> Dict[str, float]:
        kon_scale = 10.0 ** vec[0]
        koff_scale = 10.0 ** vec[1]
        pd1_50_density = 10.0 ** vec[2]
        internalization = 10.0 ** vec[3]
        params = dict(base_params)
        params["kon_PD1_PDL1"] = base_params["kon_PD1_PDL1"] * kon_scale
        params["kon_PD1_PDL2"] = base_params["kon_PD1_PDL2"] * kon_scale
        params["kon_PD1_aPD1"] = base_params["kon_PD1_aPD1"] * kon_scale
        params["koff_PD1_PDL1"] = base_params["koff_PD1_PDL1"] * koff_scale
        params["koff_PD1_PDL2"] = base_params["koff_PD1_PDL2"] * koff_scale
        params["koff_PD1_aPD1"] = base_params["koff_PD1_aPD1"] * koff_scale
        params["pd1_whitebox_pd1_50_density"] = pd1_50_density
        params["pd1_occ_internalization_per_day"] = internalization
        return params

    def objective(vec: np.ndarray) -> float:
        params = make_params(vec)
        accum = 0.0
        for scenario in scenario_list:
            preds = simulate_scenario(scenario, params)
            if not np.all(np.isfinite(preds)):
                return 1e6
            diff = preds - scenario.target
            accum += float(np.dot(diff, diff))
        return accum / max(total_points, 1)

    x0 = np.log10([1.0, 1.0, base_params.get("PD1_50", 1.0), base_params.get("pd1_occ_internalization_per_day", 0.01)])
    bounds = [
        (-4, 4),   # kon scale
        (-4, 4),   # koff scale
        (-2, 4),   # PD1_50 density
        (-6, 2),   # internalisation per day
    ]

    result = minimize(
        objective,
        x0=np.array(x0, dtype=float),
        bounds=bounds,
        method="L-BFGS-B",
    )
    best_params = make_params(result.x)
    payload = {
        "loss": float(result.fun),
        "success": bool(result.success),
        "message": result.message,
        "log10_params": result.x.tolist(),
        "scaled_parameters": {
            "kon_PD1_PDL1": best_params["kon_PD1_PDL1"],
            "kon_PD1_PDL2": best_params["kon_PD1_PDL2"],
            "kon_PD1_aPD1": best_params["kon_PD1_aPD1"],
            "koff_PD1_PDL1": best_params["koff_PD1_PDL1"],
            "koff_PD1_PDL2": best_params["koff_PD1_PDL2"],
            "koff_PD1_aPD1": best_params["koff_PD1_aPD1"],
            "pd1_whitebox_pd1_50_density": best_params["pd1_whitebox_pd1_50_density"],
            "pd1_occ_internalization_per_day": best_params["pd1_occ_internalization_per_day"],
        },
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(payload, indent=2))
    return best_params


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Fit PD-1 white-box parameters using MATLAB training data.")
    parser.add_argument("--training-path", type=Path, default=Path("artifacts/training/pd1_whitebox_training.parquet"))
    parser.add_argument("--parameter-file", type=Path, default=Path("parameters/example1_parameters.json"))
    parser.add_argument("--output-json", type=Path, default=Path("artifacts/training/pd1_whitebox_fit.json"))
    args = parser.parse_args(argv)

    if not args.training_path.is_file():
        raise SystemExit(f"Training dataset {args.training_path} is missing.")

    frame = pd.read_parquet(args.training_path)
    scenarios = prepare_training(frame)
    base_params = load_parameters(args.parameter_file)
    valid_scenarios = filter_valid_scenarios(scenarios, base_params)
    best_params = fit_parameters(valid_scenarios, base_params, output_json=args.output_json)

    print("Fit complete. Key parameters:")
    for key in [
        "kon_PD1_PDL1",
        "kon_PD1_PDL2",
        "kon_PD1_aPD1",
        "koff_PD1_PDL1",
        "koff_PD1_PDL2",
        "koff_PD1_aPD1",
        "pd1_whitebox_pd1_50_density",
        "pd1_occ_internalization_per_day",
    ]:
        print(f"  {key}: {best_params[key]:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
