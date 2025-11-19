"""White-box helper for T-cell proliferation/trafficking."""

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
class TCellParams:
    k_act: float
    k_pro: float
    k_death: float
    k_pd1: float
    k_n_death: float
    k_treg: float
    k_t0_death: float
    k_cell_clear: float
    thymus_rate: float
    q_ln_in: float
    q_ln_out: float
    q_peripheral_in: float
    q_peripheral_out: float
    q_tumour_in: float
    n_clones: float
    n_production_rate: float
    n_production_half: float
    aT_multiplier: float
    treg_relax_rate: float
    texh_relax_rate: float
    density_tau: float
    density_scale: float
    solver_rtol: float
    solver_atol: float
    max_step_days: float

    @staticmethod
    def from_mapping(parameters: Mapping[str, float]) -> "TCellParams":
        def param(key: str, default: float) -> float:
            return _finite(parameters.get(key), default)

        div_t1 = max(param("div_T1", 1.0), 1e-9)
        n_half = param("K_nT1_pro", 1e9) / div_t1
        return TCellParams(
            k_act=param("k_T1_act", 1.0),
            k_pro=param("k_T1_pro", 0.5),
            k_death=param("k_T1_death", 0.5),
            k_pd1=param("k_T1", 0.1),
            k_n_death=param("k_nT1_death", 0.01),
            k_treg=param("k_Treg", 0.0),
            k_t0_death=param("k_T0_death", 0.0),
            k_cell_clear=param("k_cell_clear", 0.0),
            thymus_rate=param("Q_nT1_thym", 1e6) / div_t1,
            q_ln_in=param("q_LN_T1_in", 0.0),
            q_ln_out=param("q_T1_LN_out", 0.0),
            q_peripheral_in=param("q_T1_P_in", 0.0),
            q_peripheral_out=param("q_T1_P_out", 0.0),
            q_tumour_in=param("q_T1_T_in", 0.0),
            n_clones=max(param("n_T1_clones", 1.0), 1.0),
            n_production_rate=param("k_nT1_pro", 0.0) / div_t1,
            n_production_half=max(n_half, 1e-6),
            aT_multiplier=2.0 ** param("N_aT", 0.0),
            treg_relax_rate=param("tcell_treg_relax_per_day", 0.1),
            texh_relax_rate=param("tcell_texh_relax_per_day", 0.1),
            density_tau=max(param("tcell_whitebox_tau_days", 2.0), 0.0),
            density_scale=max(param("tcell_density_scale", 1.0), 0.0),
            solver_rtol=max(param("tcell_whitebox_rtol", 1e-6), 1e-9),
            solver_atol=max(param("tcell_whitebox_atol", 1e-3), 1e-12),
            max_step_days=max(param("tcell_whitebox_max_step_days", 0.5), 1e-6),
        )


@dataclass
class TCellWhiteboxModel:
    params: TCellParams
    solver_config: SolverConfig
    n_central: float
    n_peripheral: float
    n_ln: float
    a_ln: float
    t_central: float
    t_peripheral: float
    t_ln: float
    t_eff: float
    t_reg: float
    t_exh: float
    density_state: float
    time_days: float
    cell_baseline: float
    v_c_l: float
    v_p_l: float
    v_ln_l: float
    treg_observable: float
    debug_last_rhs: Optional[np.ndarray] = None
    debug_last_fluxes: Optional[Mapping[str, float]] = None

    @classmethod
    def from_context(
        cls,
        parameters: Mapping[str, float],
        context: Mapping[str, float],
        *,
        solver_config: Optional[SolverConfig] = None,
    ) -> "TCellWhiteboxModel":
        params = TCellParams.from_mapping(parameters)

        def ctx(key: str, fallback: float = 0.0) -> float:
            return _finite(context.get(key), fallback)

        solver = solver_config or SolverConfig(
            method="BDF",
            rtol=params.solver_rtol,
            atol=params.solver_atol,
            max_step=params.max_step_days,
            seed=None,
        )
        return cls(
            params=params,
            solver_config=solver,
            n_central=max(ctx("V_C.nT1", ctx("V_P.nT1", 1e5)), 0.0),
            n_peripheral=max(ctx("V_P.nT1", 1e5), 0.0),
            n_ln=max(ctx("V_LN.nT1", 1e6), 0.0),
            a_ln=max(ctx("V_LN.aT1", 1e3), 0.0),
            t_central=max(ctx("V_C.T1", ctx("V_T.T1", 0.0)), 0.0),
            t_peripheral=max(ctx("V_P.T1", ctx("V_LN.T1", 0.0)), 0.0),
            t_ln=max(ctx("V_LN.T1", ctx("V_T.T1", 0.0)), 0.0),
            t_eff=max(ctx("V_T.T1", 0.0), 0.0),
            t_reg=max(ctx("V_T.T0", 0.0), 0.0),
            t_exh=max(ctx("V_T.T_exh", 0.0), 0.0),
            density_state=max(ctx("tcell_density_per_ul", 0.0), 0.0),
            time_days=ctx("time_days", ctx("t", 0.0)),
            cell_baseline=max(ctx("cell", 1.0), 1.0),
            v_c_l=max(ctx("V_C", 5.0), 1e-9),
            v_p_l=max(ctx("V_P", 1.0), 1e-9),
            v_ln_l=max(ctx("V_LN", 1.0), 1e-9),
            treg_observable=max(ctx("V_T.T0", 0.0), 0.0),
        )

    def _volume_l(self, context: Mapping[str, float]) -> float:
        vt_l = _finite(context.get("tumour_volume_l"), 0.0)
        if vt_l > 0.0:
            return vt_l
        vt_ul = _finite(context.get("V_T"), 0.0)
        if vt_ul > 0.0:
            return vt_ul * 1e-6
        return 1e-6

    def step(self, context: MutableMapping[str, float], delta_t: float, occupancy: float) -> None:
        dt = max(float(delta_t), 0.0)
        if dt <= 0.0:
            self._writeback(context, self._volume_l(context))
            return

        start_time = self.time_days
        end_time = start_time + dt
        self.time_days = end_time

        H_mAPC = _finite(context.get("H_mAPC"), 1.0)
        H_P1 = _finite(context.get("H_P1"), 1.0)
        target_treg = max(_finite(context.get("V_T.T0"), self.t_reg), 0.0)
        self.treg_observable = target_treg
        volume_l = self._volume_l(context)
        occ = min(max(float(occupancy), 0.0), 1.0)
        c_total = max(_finite(context.get("C_total"), context.get("C1")), 0.0)
        vt_l = max(volume_l, 1e-12)

        y0 = np.array(
            [
                self.n_central,
                self.n_peripheral,
                self.n_ln,
                self.a_ln,
                self.t_central,
                self.t_peripheral,
                self.t_ln,
                self.t_eff,
                self.t_exh,
            ],
            dtype=float,
        )

        def compute_rhs(state: np.ndarray) -> tuple[np.ndarray, Mapping[str, float]]:
            (
                n_c,
                n_p,
                n_ln,
                a_ln,
                t_c,
                t_p,
                t_ln,
                t_eff,
                t_exh,
            ) = state
            k_act_eff = self.params.k_act * H_mAPC * H_P1
            thymus = self.params.thymus_rate
            denom_p = self.params.n_production_half + max(n_p, 0.0)
            denom_c = self.params.n_production_half + max(n_c, 0.0)
            prod_p = 0.0 if denom_p <= 0.0 else self.params.n_production_rate * n_p / denom_p
            prod_c = 0.0 if denom_c <= 0.0 else self.params.n_production_rate * n_c / denom_c
            r5 = thymus
            r6 = prod_p
            r7 = prod_c
            r8 = self.params.k_n_death * n_p
            r9 = self.params.k_n_death * n_c
            r10 = self.params.q_peripheral_in * n_c
            r11 = self.params.q_peripheral_out * n_p
            r12 = self.params.q_ln_in * self.v_c_l * n_c
            r13 = self.params.q_ln_out * n_ln
            r14 = k_act_eff * n_ln
            r15 = k_act_eff * n_ln * self.params.n_clones
            r16 = self.params.k_pro * a_ln
            r17 = self.params.k_pro * self.params.aT_multiplier * a_ln
            r18 = self.params.k_death * t_c
            r19 = self.params.k_death * t_p
            r20 = self.params.k_death * t_eff
            r21 = self.params.k_death * t_ln
            denom_total = max(c_total + t_eff + self.t_reg + self.cell_baseline, 1e-9)
            r22 = self.params.k_treg * t_eff * self.t_reg / denom_total
            r23 = self.params.k_pd1 * t_eff * c_total / denom_total * occ
            r24 = self.params.q_peripheral_in * t_c
            r25 = self.params.q_peripheral_out * t_p
            r26 = self.params.q_tumour_in * vt_l * t_c
            r27 = self.params.q_ln_out * t_ln
            d_n_c = r5 + r7 - r9 - r10 + r11 - r12 + r13
            d_n_p = r6 - r8 + r10 - r11
            d_n_ln = r12 - r13 - r14
            d_a_ln = r15 - r16
            d_t_c = -r18 - r24 + r25 - r26 + r27
            d_t_p = -r19 + r24 - r25
            d_t_ln = r17 - r21 - r27
            d_t_eff = r26 - r20 - r22 - r23
            d_t_exh = -self.params.k_cell_clear * t_exh + r20 + r22 + r23 + self.params.k_t0_death * self.t_reg
            return (
                np.array([d_n_c, d_n_p, d_n_ln, d_a_ln, d_t_c, d_t_p, d_t_ln, d_t_eff, d_t_exh], dtype=float),
                {
                    "r5": r5,
                    "r6": r6,
                    "r7": r7,
                    "r8": r8,
                    "r9": r9,
                    "r10": r10,
                    "r11": r11,
                    "r12": r12,
                    "r13": r13,
                    "r14": r14,
                    "r15": r15,
                    "r16": r16,
                    "r17": r17,
                    "r18": r18,
                    "r19": r19,
                    "r20": r20,
                    "r21": r21,
                    "r22": r22,
                    "r23": r23,
                    "r24": r24,
                    "r25": r25,
                    "r26": r26,
                    "r27": r27,
                    "d_t_eff": d_t_eff,
                },
            )

        def rhs(_, y: np.ndarray) -> np.ndarray:
            deriv, _ = compute_rhs(y)
            return deriv

        initial_rhs, initial_flux = compute_rhs(y0)
        self.debug_last_rhs = initial_rhs
        self.debug_last_fluxes = initial_flux

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

        self.n_central = max(final[0], 0.0)
        self.n_peripheral = max(final[1], 0.0)
        self.n_ln = max(final[2], 0.0)
        self.a_ln = max(final[3], 0.0)
        self.t_central = max(final[4], 0.0)
        self.t_peripheral = max(final[5], 0.0)
        self.t_ln = max(final[6], 0.0)
        self.t_eff = max(final[7], 0.0)
        self.t_exh = max(final[8], 0.0)
        if self.params.treg_relax_rate > 0.0 and dt > 0.0:
            d_reg = self.params.treg_relax_rate * (target_treg - self.t_reg)
            self.t_reg = max(self.t_reg + d_reg * dt, 0.0)
        else:
            self.t_reg = max(target_treg, 0.0)

        self._update_density(volume_l, dt)
        self._writeback(context, volume_l)

    def _update_density(self, volume_l: float, delta_t: float) -> None:
        denom = max(volume_l * 1e6, 1e-9)
        total_tcells = max(self.t_eff + self.treg_observable, 0.0)
        density_raw = max(total_tcells * self.params.density_scale / denom, 0.0)
        if self.params.density_tau > 0.0:
            alpha = math.exp(-delta_t / self.params.density_tau)
            self.density_state = density_raw + (self.density_state - density_raw) * alpha
        else:
            self.density_state = density_raw

    def _writeback(self, context: MutableMapping[str, float], volume_l: float) -> None:
        context["V_C.nT1"] = self.n_central
        context["V_P.nT1"] = self.n_peripheral
        context["V_LN.nT1"] = self.n_ln
        context["V_LN.aT1"] = self.a_ln
        context["V_C.T1"] = self.t_central
        context["V_P.T1"] = self.t_peripheral
        context["V_LN.T1"] = self.t_ln
        context["V_T.T1"] = self.t_eff
        context["V_T.T0"] = self.t_reg
        context["V_T.T_exh"] = self.t_exh
        context["tcell_alignment_state"] = self.t_eff
        context["tcell_density_per_ul"] = max(self.density_state, 0.0)
        context["tcell_kill_rate_hat"] = float(context.get("k_C_T1", 0.0)) * self.t_eff
