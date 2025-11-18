#!/usr/bin/env python
"""Dev helper to verify PD-1 synapse steady-state alignment.

Modes:
  * pdl1        – only Reaction 89 (PD1 + PDL1).
  * pdl1_pdl2   – add Reaction 90 (PDL2).
  * with_ab     – add Reaction 91 (aPD1 binding, no dimer).
  * with_dimer  – add Reaction 92 (cross-arm dimer).

For each mode we solve the algebraic steady state analytically, then
integrate the PD1WhiteboxModel to the same configuration.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import replace
from pathlib import Path
from typing import Iterable, Tuple

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.offline.modules.pd1_params import PD1Params, load_pd1_parameters_from_file
from src.offline.modules.pd1_whitebox import PD1WhiteboxModel
from src.offline.segment_integrator import SolverConfig

SECONDS_PER_DAY = 86400.0


def _trim_params(
    params: PD1Params,
    *,
    include_pdl2: bool,
    include_ab: bool,
    include_dimer: bool,
) -> PD1Params:
    return replace(
        params,
        kon_pd1_pdl2=params.kon_pd1_pdl2 if include_pdl2 else 0.0,
        koff_pd1_pdl2=params.koff_pd1_pdl2 if include_pdl2 else 0.0,
        kon_pd1_ab=params.kon_pd1_ab if include_ab else 0.0,
        koff_pd1_ab=params.koff_pd1_ab if include_ab else 0.0,
        chi_pd1=params.chi_pd1 if include_dimer else 0.0,
        internalisation_per_day=0.0,
        total_pdl2_molecules=params.total_pdl2_molecules if include_pdl2 else 0.0,
    )


def _c_from_alpha(alpha: float, pd1_free: float, ligand_total: float) -> float:
    if alpha <= 0.0 or pd1_free <= 0.0 or ligand_total <= 0.0:
        return 0.0
    denom = 1.0 + alpha * pd1_free
    if denom <= 0.0:
        return 0.0
    return (alpha * pd1_free * ligand_total) / denom


def _ab_density(
    pd1_free: float,
    params: PD1Params,
    ab_effective: float,
    include_dimer: bool,
) -> float:
    if params.kon_pd1_ab <= 0.0 or pd1_free <= 0.0 or ab_effective <= 0.0:
        return 0.0
    denom = params.koff_pd1_ab + params.internalisation_per_day
    if include_dimer:
        denom += params.chi_pd1 * params.kon_pd1_ab * pd1_free
    if denom <= 0.0:
        return 0.0
    value = (2.0 * params.kon_pd1_ab * pd1_free * ab_effective) / denom
    return max(value, 0.0)


def _dimer_density(pd1_free: float, ab_density: float, params: PD1Params) -> float:
    if pd1_free <= 0.0 or ab_density <= 0.0 or params.kon_pd1_ab <= 0.0:
        return 0.0
    denom = 2.0 * params.koff_pd1_ab + params.internalisation_per_day
    if denom <= 0.0:
        return 0.0
    rate = params.chi_pd1 * params.kon_pd1_ab * pd1_free * ab_density
    return max(rate / denom, 0.0)


def _solve_pd1_free(
    params: PD1Params,
    pd1_total: float,
    *,
    alpha_pdl1: float,
    pdl1_total: float,
    alpha_pdl2: float,
    pdl2_total: float,
    ab_effective: float,
    include_pdl2: bool,
    include_ab: bool,
    include_dimer: bool,
) -> float:
    def total_from_free(x: float) -> float:
        total = x
        total += _c_from_alpha(alpha_pdl1, x, pdl1_total)
        if include_pdl2:
            total += _c_from_alpha(alpha_pdl2, x, pdl2_total)
        if include_ab:
            ab_density = _ab_density(x, params, ab_effective, include_dimer)
            total += ab_density
            if include_dimer:
                dimer_density = _dimer_density(x, ab_density, params)
                total += 2.0 * dimer_density
        return total

    lo = 0.0
    hi = max(pd1_total, 1e-12)
    target = pd1_total
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        err = total_from_free(mid) - target
        if abs(err) <= 1e-12:
            return mid
        if err > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def analytic_state(
    params: PD1Params,
    *,
    pd1_density: float,
    pdl1_density: float,
    pdl2_density: float,
    ab_molar: float,
    include_pdl2: bool,
    include_ab: bool,
    include_dimer: bool,
) -> Tuple[float, float, float, float]:
    alpha_pdl1 = params.kon_pd1_pdl1 / max(params.koff_pd1_pdl1, 1e-18)
    alpha_pdl2 = params.kon_pd1_pdl2 / max(params.koff_pd1_pdl2, 1e-18) if include_pdl2 else 0.0
    ab_effective = ab_molar / max(params.gamma_c_nivolumab, 1e-18)
    pd1_free = _solve_pd1_free(
        params,
        pd1_density,
        alpha_pdl1=alpha_pdl1,
        pdl1_total=pdl1_density,
        alpha_pdl2=alpha_pdl2 if include_pdl2 else 0.0,
        pdl2_total=pdl2_density if include_pdl2 else 0.0,
        ab_effective=ab_effective,
        include_pdl2=include_pdl2,
        include_ab=include_ab,
        include_dimer=include_dimer,
    )
    c1 = _c_from_alpha(alpha_pdl1, pd1_free, pdl1_density)
    c2 = _c_from_alpha(alpha_pdl2, pd1_free, pdl2_density) if include_pdl2 else 0.0
    ab_complex = _ab_density(pd1_free, params, ab_effective, include_dimer) if include_ab else 0.0
    dimer = _dimer_density(pd1_free, ab_complex, params) if include_ab and include_dimer else 0.0
    return c1, c2, ab_complex, dimer


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Verify PD1 synapse steady state between analytic solution and ODE.")
    parser.add_argument(
        "--parameter-file",
        type=Path,
        default=Path("parameters/example1_parameters.json"),
        help="Path to parameter JSON (default: parameters/example1_parameters.json)",
    )
    parser.add_argument(
        "--mode",
        choices=["pdl1", "pdl1_pdl2", "with_ab", "with_dimer"],
        default="pdl1",
        help="System configuration to evaluate (default: pdl1).",
    )
    parser.add_argument(
        "--pd1-density",
        type=float,
        help="Total PD1 density (molecules/µm^2). Defaults to PD1Params.pd1_surface_density().",
    )
    parser.add_argument(
        "--pdl1-density",
        type=float,
        help="Total PDL1 density (molecules/µm^2). Defaults to PD1Params.pdl1_surface_density().",
    )
    parser.add_argument(
        "--pdl2-density",
        type=float,
        help="Total PDL2 density (molecules/µm^2). Defaults to PD1Params.pdl2_surface_density().",
    )
    parser.add_argument(
        "--kon",
        type=float,
        help="Override kon for PD1-PDL1 (1/(density*day)). Defaults to PD1Params.kon_pd1_pdl1.",
    )
    parser.add_argument(
        "--koff",
        type=float,
        help="Override koff for PD1-PDL1 (1/day). Defaults to PD1Params.koff_pd1_pdl1.",
    )
    parser.add_argument(
        "--sim-seconds",
        type=float,
        default=1e4,
        help="Total integration time in seconds (default: 1e4). Converted internally to days.",
    )
    parser.add_argument(
        "--chunk-seconds",
        type=float,
        default=1000.0,
        help="Integration chunk size in seconds (default: 1000).",
    )
    parser.add_argument(
        "--ab-molar",
        type=float,
        default=0.0,
        help="Constant aPD1 concentration (molar) to feed into the white-box.",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        help="Override solver relative tolerance.",
    )
    parser.add_argument(
        "--atol",
        type=float,
        help="Override solver absolute tolerance.",
    )
    return parser.parse_args(argv)


def integrate_to_equilibrium(
    params: PD1Params,
    *,
    pd1_density: float,
    pdl1_density: float,
    pdl2_density: float,
    total_days: float,
    chunk_days: float,
    solver_cfg: SolverConfig,
    ab_molar: float,
) -> PD1WhiteboxModel:
    context = {
        "syn_pd1_total": pd1_density,
        "syn_pdl1_total": pdl1_density,
        "syn_pdl2_total": pdl2_density,
        "syn_pd1_pdl1": 0.0,
        "syn_pd1_pdl2": 0.0,
        "syn_pd1_apd1": 0.0,
        "syn_pd1_apd1_pd1": 0.0,
        "time_days": 0.0,
    }
    model = PD1WhiteboxModel.from_context(params, context, solver_config=solver_cfg)
    remaining = max(total_days, 0.0)
    while remaining > 0.0:
        dt = min(remaining, chunk_days if chunk_days > 0 else remaining)
        model.step(ab_molar, dt)
        remaining -= dt
    return model


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    params = load_pd1_parameters_from_file(args.parameter_file)
    pd1_density = float(args.pd1_density if args.pd1_density is not None else params.pd1_surface_density())
    pdl1_density = float(args.pdl1_density if args.pdl1_density is not None else params.pdl1_surface_density())
    pdl2_density = float(args.pdl2_density if args.pdl2_density is not None else params.pdl2_surface_density())
    kon = float(args.kon if args.kon is not None else params.kon_pd1_pdl1)
    koff = float(args.koff if args.koff is not None else params.koff_pd1_pdl1)
    include_pdl2 = args.mode in ("pdl1_pdl2", "with_ab", "with_dimer")
    include_ab = args.mode in ("with_ab", "with_dimer")
    include_dimer = args.mode == "with_dimer"

    trimmed = _trim_params(
        replace(params, kon_pd1_pdl1=kon, koff_pd1_pdl1=koff),
        include_pdl2=include_pdl2,
        include_ab=include_ab,
        include_dimer=include_dimer,
    )

    total_days = max(args.sim_seconds, 0.0) / SECONDS_PER_DAY
    chunk_days = max(args.chunk_seconds, 1.0) / SECONDS_PER_DAY
    solver_cfg = SolverConfig(
        method="BDF",
        rtol=float(args.rtol if args.rtol is not None else trimmed.solver_rtol),
        atol=float(args.atol if args.atol is not None else trimmed.solver_atol),
        max_step=trimmed.max_step_days,
        seed=None,
    )

    numeric_model = integrate_to_equilibrium(
        trimmed,
        pd1_density=pd1_density,
        pdl1_density=pdl1_density,
        pdl2_density=pdl2_density if include_pdl2 else 0.0,
        total_days=total_days,
        chunk_days=chunk_days,
        solver_cfg=solver_cfg,
        ab_molar=float(args.ab_molar if include_ab else 0.0),
    )

    analytic_pdl1, analytic_pdl2, analytic_ab, analytic_dimer = analytic_state(
        trimmed,
        pd1_density=pd1_density,
        pdl1_density=pdl1_density,
        pdl2_density=pdl2_density if include_pdl2 else 0.0,
        ab_molar=float(args.ab_molar if include_ab else 0.0),
        include_pdl2=include_pdl2,
        include_ab=include_ab,
        include_dimer=include_dimer,
    )

    def _report(name: str, numeric: float, analytic: float) -> None:
        diff = numeric - analytic
        rel = diff / analytic if analytic not in (0.0, -0.0) else float("nan")
        print(f"{name:>15}: numeric={numeric: .9g} analytic={analytic: .9g} diff={diff:+.3e} rel={rel:+.3e}")

    print("=== PD1 steady-state check ===")
    print(f"mode: {args.mode}")
    print(f"pd1_density: {pd1_density:.6g} molecules/µm^2")
    print(f"pdl1_density: {pdl1_density:.6g} molecules/µm^2")
    if include_pdl2:
        print(f"pdl2_density: {pdl2_density:.6g} molecules/µm^2")
    if include_ab:
        print(f"aPD1_molar: {args.ab_molar:.6g} M")
    print(f"kon: {kon:.6g} 1/(density*day)")
    print(f"koff: {koff:.6g} 1/day")
    _report("PD1_PDL1", numeric_model.syn_pd1_pdl1, analytic_pdl1)
    if include_pdl2:
        _report("PD1_PDL2", numeric_model.syn_pd1_pdl2, analytic_pdl2)
    if include_ab:
        _report("PD1_aPD1", numeric_model.syn_pd1_ab, analytic_ab)
    if include_dimer:
        _report("PD1_aPD1_PD1", numeric_model.syn_pd1_ab_pd1, analytic_dimer)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
