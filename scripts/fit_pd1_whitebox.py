#!/usr/bin/env python

"""Fit PD-1 white-box parameters against MATLAB-generated training data."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.offline.modules.pd1_params import PD1Params, load_pd1_parameters_from_file
from src.offline.segment_integrator import SolverConfig
from src.offline.modules.pd1_whitebox import PD1WhiteboxModel


@dataclass
class ScenarioData:
    name: str
    dt: np.ndarray
    concentration: np.ndarray
    target: np.ndarray
    init_context: Dict[str, float]


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
    parameters: PD1Params,
) -> np.ndarray:
    context = dict(scenario.init_context)
    solver_cfg = SolverConfig(method="BDF", rtol=parameters.solver_rtol, atol=parameters.solver_atol, max_step=parameters.max_step_days, seed=None)
    model = PD1WhiteboxModel.from_context(parameters, context, solver_config=solver_cfg)
    preds = np.empty_like(scenario.target)
    for idx, (dt, conc) in enumerate(zip(scenario.dt, scenario.concentration)):
        outputs = model.step(conc, max(dt, 0.0))
        preds[idx] = outputs.occupancy
    return preds


def filter_valid_scenarios(
    scenarios: Iterable[ScenarioData],
    parameters: PD1Params,
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
    base_params: PD1Params,
    *,
    output_json: Path,
) -> PD1Params:
    scenario_list = list(scenarios)
    total_points = sum(len(s.target) for s in scenario_list)

    def make_params(vec: np.ndarray) -> PD1Params:
        kon_scale = 10.0 ** vec[0]
        koff_scale = 10.0 ** vec[1]
        pd1_50_density = 10.0 ** vec[2]
        internalization = 10.0 ** vec[3]
        return replace(
            base_params,
            kon_pd1_pdl1=base_params.kon_pd1_pdl1 * kon_scale,
            kon_pd1_pdl2=base_params.kon_pd1_pdl2 * kon_scale,
            kon_pd1_ab=base_params.kon_pd1_ab * kon_scale,
            koff_pd1_pdl1=base_params.koff_pd1_pdl1 * koff_scale,
            koff_pd1_pdl2=base_params.koff_pd1_pdl2 * koff_scale,
            koff_pd1_ab=base_params.koff_pd1_ab * koff_scale,
            pd1_50_density=pd1_50_density,
            internalisation_per_day=internalization,
        )

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
        "scaled_parameters": asdict(best_params),
    }
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(payload, indent=2))
    return best_params


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Fit PD-1 white-box parameters using MATLAB training data.")
    parser.add_argument("--training-path", type=Path, default=Path("artifacts/training/pd1_whitebox_training.parquet"))
    parser.add_argument("--parameter-file", type=Path, default=Path("parameters/example1_parameters.json"))
    parser.add_argument("--output-json", type=Path, default=Path("artifacts/training/pd1_whitebox_fit.json"))
    parser.add_argument(
        "--max-scenarios",
        type=int,
        default=None,
        help="Optional cap on the number of training scenarios to accelerate local experiments.",
    )
    parser.add_argument(
        "--pd1-max-step-days",
        type=float,
        default=0.05,
        help="Override pd1_whitebox_max_step_days during fitting (default: 0.05 days).",
    )
    args = parser.parse_args(argv)

    if not args.training_path.is_file():
        raise SystemExit(f"Training dataset {args.training_path} is missing.")

    frame = pd.read_parquet(args.training_path)
    scenarios = prepare_training(frame)
    if args.max_scenarios is not None and args.max_scenarios > 0:
        scenarios = scenarios[: args.max_scenarios]
    base_params = load_pd1_parameters_from_file(args.parameter_file)
    if args.pd1_max_step_days is not None and args.pd1_max_step_days > 0:
        base_params = replace(base_params, max_step_days=float(args.pd1_max_step_days))
    valid_scenarios = filter_valid_scenarios(scenarios, base_params)
    best_params = fit_parameters(valid_scenarios, base_params, output_json=args.output_json)

    print("Fit complete. Key parameters:")
    key_map = {
        "kon_PD1_PDL1": best_params.kon_pd1_pdl1,
        "kon_PD1_PDL2": best_params.kon_pd1_pdl2,
        "kon_PD1_aPD1": best_params.kon_pd1_ab,
        "koff_PD1_PDL1": best_params.koff_pd1_pdl1,
        "koff_PD1_PDL2": best_params.koff_pd1_pdl2,
        "koff_PD1_aPD1": best_params.koff_pd1_ab,
        "pd1_whitebox_pd1_50_density": best_params.pd1_50_density,
        "pd1_occ_internalization_per_day": best_params.internalisation_per_day,
    }
    for key, value in key_map.items():
        print(f"  {key}: {value:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
