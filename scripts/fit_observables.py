"""Grey-box fitting utility for snapshot observables.

This script tunes a small set of runtime parameters (e.g., pd1 occupancy filter
and dynamic geometry constants) so that the Python surrogate matches the
MATLAB reference trajectories for a chosen scenario.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from scripts.validate_surrogate import SCENARIO_REGISTRY
from src.offline.entities import ScenarioResult
from src.offline.simulation import simulate_frozen_model

DEFAULT_MODULE_BLOCKS: Sequence[str] = (
    "pd1_bridge_block",
    "pd1_occupancy_filter_block",
    "tumour_geometry_dynamic_block",
)
DEFAULT_OBSERVABLES: Sequence[str] = (
    "tumour_volume_l",
    "pd1_occupancy",
    "tcell_density_per_ul",
)
ARTIFACTS_DIR = Path("artifacts/validation")


@dataclass
class ParamSpec:
    name: str
    initial: float
    lower: float
    upper: float


def _parse_param_spec(arg: str) -> ParamSpec:
    """Parse --param name=init:lower:upper."""

    if "=" not in arg:
        raise argparse.ArgumentTypeError("expected name=value[:lower[:upper]]")
    name, payload = arg.split("=", 1)
    parts = payload.split(":")
    try:
        initial = float(parts[0])
    except ValueError as exc:  # pragma: no cover - CLI validation
        raise argparse.ArgumentTypeError(f"invalid initial value in '{arg}': {exc}") from exc
    lower = float(parts[1]) if len(parts) > 1 and parts[1] else -math.inf
    upper = float(parts[2]) if len(parts) > 2 and parts[2] else math.inf
    if lower > initial or initial > upper:
        raise argparse.ArgumentTypeError(
            f"initial value {initial} not within bounds [{lower}, {upper}]"
        )
    return ParamSpec(name=name.strip(), initial=initial, lower=lower, upper=upper)


def _load_reference_frame(scenario: str) -> pd.DataFrame:
    path = ARTIFACTS_DIR / f"{scenario}_reference.csv"
    if not path.is_file():
        raise FileNotFoundError(
            f"Reference trajectory {path} not found; generate it via MATLAB exporter"
        )
    return pd.read_csv(path, comment="#")


def _load_override_json(path: Path) -> Dict[str, float]:
    try:
        payload = json.loads(path.read_text(encoding="utf8"))
    except Exception as exc:  # pragma: no cover - CLI validation
        raise SystemExit(f"Failed to load overrides from '{path}': {exc}") from exc
    results: Dict[str, float] = {}
    if isinstance(payload, dict):
        items = payload.items()
    elif isinstance(payload, list):
        items = []
        for entry in payload:
            if isinstance(entry, dict) and "name" in entry and "value" in entry:
                items.append((entry["name"], entry["value"]))
            else:
                raise SystemExit(
                    f"Invalid override entry in '{path}': {entry!r}; expected objects with 'name' and 'value'"
                )
    else:
        raise SystemExit(f"Invalid override payload in '{path}'; expected dict or list")
    for name, value in items:
        try:
            results[str(name).strip()] = float(value)
        except ValueError as exc:
            raise SystemExit(
                f"Non-numeric override value for '{name}' in '{path}': {exc}"
            ) from exc
    return results


def _simulate_for_params(
    scenario_name: str,
    overrides: Dict[str, float],
    module_blocks: Sequence[str],
    ic_mode: str,
    ic_options,
) -> ScenarioResult:
    scenario = SCENARIO_REGISTRY[scenario_name]
    simulate_kwargs = {
        "sample_interval_hours": scenario.sample_interval_hours,
        "custom_doses": scenario.custom_doses,
        "context_outputs": scenario.context_outputs,
        "module_blocks": module_blocks,
        "param_overrides": overrides,
        "ic_mode": ic_mode,
    }
    if ic_mode == "target_volume":
        simulate_kwargs["ic_options"] = ic_options
    return simulate_frozen_model(
        scenario.snapshot,
        days=scenario.stop_time,
        therapy=scenario.therapy,
        emit_diagnostics=False,
        **simulate_kwargs,
    )


def _merge_frames(result: ScenarioResult, reference: pd.DataFrame) -> pd.DataFrame:
    surrogate = result.to_frame()
    return surrogate.merge(reference, on="time_days", suffixes=("_sur", "_ref"))


def _build_residuals(
    merged: pd.DataFrame,
    observables: Sequence[str],
    min_ref: float,
    weights: Dict[str, float],
) -> np.ndarray:
    residuals: List[np.ndarray] = []
    for column in observables:
        sur = merged[f"{column}_sur"].to_numpy(dtype=float)
        ref = merged[f"{column}_ref"].to_numpy(dtype=float)
        denom = np.maximum(np.abs(ref), min_ref)
        diff = (sur - ref) / denom
        weight = math.sqrt(weights.get(column, 1.0))
        residuals.append(weight * diff)
    return np.concatenate(residuals, axis=0)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Fit grey-box observables parameters")
    parser.add_argument("--scenario", default="A1", choices=sorted(SCENARIO_REGISTRY.keys()))
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="name=init[:lower[:upper]]",
        help="Parameter to fit (may repeat)",
        dest="params",
    )
    parser.add_argument(
        "--observable",
        action="append",
        default=list(DEFAULT_OBSERVABLES),
        help="Observable to include (default: tumour_volume_l, pd1_occupancy, tcell_density_per_ul)",
    )
    parser.add_argument(
        "--observable-weight",
        action="append",
        default=[],
        metavar="name=weight",
        help="Optional weight per observable",
    )
    parser.add_argument(
        "--max-iter", type=int, default=8, help="Maximum optimizer iterations"
    )
    parser.add_argument(
        "--min-ref-value",
        type=float,
        default=1e-6,
        help="Floor applied to reference magnitudes when computing residuals",
    )
    parser.add_argument(
        "--ic-mode",
        default="snapshot",
        choices=["snapshot", "target_volume"],
        help="Initial condition mode passed to simulate_frozen_model",
    )
    parser.add_argument(
        "--ic-target-diam-cm",
        type=float,
        default=0.5,
        help="Target diameter (cm) when ic-mode=target_volume",
    )
    parser.add_argument(
        "--ic-max-days",
        type=float,
        default=150.0,
        help="Max days for target-volume IC generation",
    )
    parser.add_argument(
        "--ic-max-wall-seconds",
        type=float,
        default=20.0,
        help="Wall-clock budget for target-volume IC generation",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        help="Optional path to write the fitted parameter overrides as JSON",
    )
    parser.add_argument(
        "--module-block",
        action="append",
        default=list(DEFAULT_MODULE_BLOCKS),
        help="Module blocks to activate during fitting",
    )
    parser.add_argument(
        "--override",
        action="append",
        default=[],
        metavar="name=value",
        help="Additional fixed param overrides (applied on top of fitted params)",
    )

    args = parser.parse_args(argv)

    if not args.params:
        parser.error("at least one --param specification is required")

    param_specs = [_parse_param_spec(p) for p in args.params]
    base_overrides: Dict[str, float] = {}
    for item in args.override:
        if not item:
            continue
        if item.startswith("@"):
            path = Path(item[1:])
            base_overrides.update(_load_override_json(path))
            continue
        if "=" not in item:
            parser.error(f"invalid --override '{item}' (expected name=value or @file)")
        name, value = item.split("=", 1)
        base_overrides[name.strip()] = float(value)

    weights: Dict[str, float] = {}
    for item in args.observable_weight:
        if "=" not in item:
            parser.error(f"invalid --observable-weight '{item}' (expected name=value)")
        name, value = item.split("=", 1)
        weights[name.strip()] = float(value)

    reference_frame = _load_reference_frame(args.scenario)
    observables = list(dict.fromkeys(args.observable))

    ic_options = None
    if args.ic_mode == "target_volume":
        from src.offline.initial_conditions import ICOptions

        ic_options = ICOptions(
            target_diameter_cm=args.ic_target_diam_cm,
            max_days=args.ic_max_days,
            max_wall_seconds=args.ic_max_wall_seconds,
            reset_policy="cancer_only",
        )

    theta0 = np.array([spec.initial for spec in param_specs], dtype=float)
    lower = np.array([spec.lower for spec in param_specs], dtype=float)
    upper = np.array([spec.upper for spec in param_specs], dtype=float)

    def residuals(theta: np.ndarray) -> np.ndarray:
        overrides = dict(base_overrides)
        overrides.update({spec.name: float(value) for spec, value in zip(param_specs, theta)})
        try:
            result = _simulate_for_params(
                args.scenario,
                overrides=overrides,
                module_blocks=args.module_block,
                ic_mode=args.ic_mode,
                ic_options=ic_options,
            )
        except Exception as exc:  # pragma: no cover - robust optimisation
            print(f"[fit] simulation failed: {exc}")
            return np.ones(len(observables) * len(reference_frame), dtype=float) * 1e3

        merged = _merge_frames(result, reference_frame)
        return _build_residuals(
            merged,
            observables=observables,
            min_ref=args.min_ref_value,
            weights=weights,
        )

    res = least_squares(
        residuals,
        theta0,
        bounds=(lower, upper),
        max_nfev=args.max_iter,
        verbose=2,
    )

    best_overrides = dict(base_overrides)
    for spec, value in zip(param_specs, res.x):
        best_overrides[spec.name] = float(value)

    print("\nBest-fit parameters:")
    for name, value in best_overrides.items():
        print(f"  {name} = {value:.6g}")
    print(f"optimizer status: {res.status} ({res.message})")
    print(f"final cost: {res.cost}")

    if args.output_json:
        args.output_json.write_text(json.dumps(best_overrides, indent=2), encoding="utf8")
        print(f"wrote overrides to {args.output_json}")

    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry
    raise SystemExit(main())
