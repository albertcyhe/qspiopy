"""White-box helper for tumour volume dynamics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, MutableMapping, Optional

import math
import numpy as np

from ..segment_integrator import SolverConfig
from ..stiff_ode import integrate_local_system


def _finite(value: object, default: float = 0.0) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


@dataclass(frozen=True)
class GeometryParams:
    k_growth: float
    k_death: float
    k_kill: float
    k_clear: float
    volume_min_l: float
    c_max: float
    vol_cell_l: float
    vol_tcell_l: float
    solver_rtol: float
    solver_atol: float
    max_step_days: float

    @staticmethod
    def from_mapping(parameters: Mapping[str, float]) -> "GeometryParams":
        def param(key: str, default: float) -> float:
            return _finite(parameters.get(key), default)

        return GeometryParams(
            k_growth=param("k_C1_growth", 0.0),
            k_death=param("k_C1_death", 0.0),
            k_kill=param("k_C_T1", 0.0),
            k_clear=param("k_cell_clear", 0.0),
            volume_min_l=max(param("V_Tmin", 0.0), 0.0),
            c_max=max(param("C_max", 1.0), 1.0),
            vol_cell_l=max(param("vol_cell", 0.0), 0.0),
            vol_tcell_l=max(param("vol_Tcell", 0.0), 0.0),
            solver_rtol=max(param("geometry_whitebox_rtol", 1e-7), 1e-12),
            solver_atol=max(param("geometry_whitebox_atol", 1e-2), 1e-12),
            max_step_days=max(param("geometry_whitebox_max_step_days", 0.25), 1e-6),
        )


@dataclass
class GeometryWhiteboxModel:
    params: GeometryParams
    solver_config: SolverConfig
    c_live: float
    c_dead: float
    volume_l: float
    time_days: float
    cell_baseline: float

    @classmethod
    def from_context(
        cls,
        parameters: Mapping[str, float],
        context: Mapping[str, float],
        *,
        solver_config: Optional[SolverConfig] = None,
    ) -> "GeometryWhiteboxModel":
        params = GeometryParams.from_mapping(parameters)

        def ctx(key: str, fallback: float = 0.0) -> float:
            return _finite(context.get(key), fallback)

        solver = solver_config or SolverConfig(
            method="BDF",
            rtol=params.solver_rtol,
            atol=params.solver_atol,
            max_step=params.max_step_days,
            seed=None,
        )
        volume = ctx("tumour_volume_l", ctx("V_T", params.volume_min_l))
        return cls(
            params=params,
            solver_config=solver,
            c_live=max(ctx("V_T.C1", ctx("C1", 0.0)), 0.0),
            c_dead=max(ctx("V_T.C_x", ctx("C_x", 0.0)), 0.0),
            volume_l=max(volume, params.volume_min_l),
            time_days=ctx("time_days", ctx("t", 0.0)),
            cell_baseline=max(ctx("cell", 1.0), 1.0),
        )

    def step(self, context: MutableMapping[str, float], delta_t: float, occupancy: float) -> None:
        dt = max(float(delta_t), 0.0)
        if dt <= 0.0:
            self._writeback(context)
            return

        start_time = self.time_days
        end_time = start_time + dt
        self.time_days = end_time

        t_eff = max(_finite(context.get("V_T.T1"), 0.0), 0.0)
        t_reg = max(_finite(context.get("V_T.T0"), 0.0), 0.0)
        total_tcells = t_eff + t_reg
        occ = min(max(float(occupancy), 0.0), 1.0)

        y0 = np.array([self.c_live, self.c_dead], dtype=float)

        def rhs(_, y: np.ndarray) -> np.ndarray:
            c_live, c_dead = y
            logistic = self.params.k_growth * c_live * max(1.0 - c_live / self.params.c_max, 0.0)
            death = self.params.k_death * c_live
            denom = max(c_live + total_tcells + self.cell_baseline, 1e-9)
            kill = self.params.k_kill * t_eff * c_live / denom * max(1.0 - occ, 0.0)
            d_live = logistic - death - kill
            d_dead = -self.params.k_clear * c_dead + death + kill
            return np.array([d_live, d_dead], dtype=float)

        try:
            sol = integrate_local_system(
                rhs,
                y0,
                start_time,
                end_time,
                solver=self.solver_config,
                max_internal_step_days=min(dt, self.params.max_step_days),
            )
            final = sol if isinstance(sol, np.ndarray) else np.asarray(sol, dtype=float)
        except Exception:
            final = y0 + rhs(0.0, y0) * dt

        self.c_live = max(final[0], 0.0)
        self.c_dead = max(final[1], 0.0)
        self._writeback(context)

    def _writeback(self, context: MutableMapping[str, float]) -> None:
        t_eff = max(_finite(context.get("V_T.T1"), 0.0), 0.0)
        t_reg = max(_finite(context.get("V_T.T0"), 0.0), 0.0)
        t_exh = max(_finite(context.get("V_T.T_exh"), context.get("T_exh", 0.0)), 0.0)
        volume = (
            self.params.volume_min_l
            + self.params.vol_cell_l * (self.c_dead + self.c_live)
            + self.params.vol_tcell_l * (t_exh + t_eff + t_reg)
        )
        self.volume_l = max(volume, self.params.volume_min_l)

        context["V_T.C1"] = self.c_live
        context["C1"] = self.c_live
        context["V_T.C_x"] = self.c_dead
        context["C_x"] = self.c_dead
        context["C_total"] = self.c_live
        context["tumour_volume_l"] = self.volume_l
        context["tumor_volume_l"] = self.volume_l
        context["V_T"] = self.volume_l * 1e6

        density = 0.0
        if self.volume_l > 0.0:
            density = max(_finite(context.get("V_T.T1"), 0.0), 0.0) / (self.volume_l * 1e6)
        context["tcell_density_per_ul"] = density

        context["geom_live_volume_l"] = self.params.vol_cell_l * self.c_live
        context["geom_dead_instant_volume_l"] = self.params.vol_cell_l * self.c_dead
        context["geom_dead_filtered_volume_l"] = context["geom_dead_instant_volume_l"]
        context["geom_tcell_volume_l"] = self.params.vol_tcell_l * (t_eff + t_reg + t_exh)
        context["geom_target_volume_l"] = self.volume_l
        context["geom_volume_smoothed_l"] = self.volume_l
        context["geom_volume_time_days"] = self.time_days


__all__ = ["GeometryWhiteboxModel", "GeometryParams"]
