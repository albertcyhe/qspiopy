"""White-box helper for nT1/aT1/T1 dynamics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, MutableMapping

import math


@dataclass
class TCellWhiteboxModel:
    """Very small ODE system for T-cell proliferation/trafficking."""

    nT1: float
    aT1: float
    t_tumour: float
    t_ln: float
    t_peripheral: float
    t_central: float
    t_treg: float

    q_ln_in: float
    q_ln_out: float
    q_peripheral_in: float
    q_peripheral_out: float
    q_tumour_in: float
    k_act: float
    k_pro: float
    k_death: float
    k_pd1_supp: float
    k_exhaustion: float
    thymus_influx: float
    carrying_capacity: float
    tcell_volume_factor: float
    smoothing_tau: float
    density_state: float

    @classmethod
    def from_context(cls, params: Mapping[str, float], context: Mapping[str, float]) -> "TCellWhiteboxModel":
        def ctx(key: str, fallback: float = 0.0) -> float:
            try:
                return float(context.get(key, fallback))
            except (TypeError, ValueError):
                return fallback

        return cls(
            nT1=ctx("tcell_ln", ctx("V_LN.nT1", 1e6)),
            aT1=ctx("tcell_ln", 0.0) * 0.01,
            t_tumour=ctx("tcell_tumour", ctx("V_T.T1", 0.0)),
            t_ln=ctx("tcell_ln", ctx("V_LN.T1", 0.0)),
            t_peripheral=ctx("tcell_peripheral", ctx("V_P.T1", 0.0)),
            t_central=ctx("tcell_central", ctx("V_C.T1", 0.0)),
            t_treg=ctx("tcell_treg", ctx("V_T.T0", 0.0)),
            q_ln_in=float(params.get("q_LN_T1_in", 0.0)),
            q_ln_out=float(params.get("q_T1_LN_out", 0.0)),
            q_peripheral_in=float(params.get("q_T1_P_in", 0.0)),
            q_peripheral_out=float(params.get("q_T1_P_out", 0.0)),
            q_tumour_in=float(params.get("q_T1_T_in", 0.0)),
            k_act=float(params.get("k_T1_act", 0.0)),
            k_pro=float(params.get("k_T1_pro", 0.0)),
            k_death=float(params.get("k_T1_death", 0.0)),
            k_pd1_supp=float(params.get("k_T1", 0.0)),
            k_exhaustion=float(params.get("k_T1", 0.0)),
            thymus_influx=float(params.get("Q_nT1_thym", 0.0)),
            carrying_capacity=max(float(params.get("T1_tumour_capacity", 8e6)), 1e5),
            tcell_volume_factor=float(params.get("geometry_tcell_volume_factor", 1.0)),
            smoothing_tau=max(float(params.get("tcell_whitebox_tau_days", 3.0)), 0.0),
            density_state=ctx("tcell_density_per_ul", 0.0),
        )

    def _logistic(self, value: float, gain: float, loss: float, dt: float) -> float:
        return max(value + (gain - loss * value) * dt, 0.0)

    def step(self, context: MutableMapping[str, float], delta_t: float, occupancy: float) -> None:
        dt = max(float(delta_t), 0.0)
        if dt <= 0.0:
            return

        suppress = max(1.0 - occupancy, 0.0)

        self.nT1 = self._logistic(self.nT1, self.thymus_influx, self.q_ln_in + self.q_ln_out + self.k_death, dt)
        activation = self.k_act * self.t_ln * suppress
        self.aT1 = self._logistic(self.aT1, activation, self.k_death, dt)

        tumour_influx = self.q_tumour_in * (self.t_peripheral + self.t_central)
        logistic = (1.0 - self.t_tumour / self.carrying_capacity)
        tumour_growth = self.k_pro * self.aT1 * max(logistic, 0.0)
        tumour_loss = (self.k_death + self.k_pd1_supp * (1.0 - suppress) + self.k_exhaustion) * self.t_tumour
        self.t_tumour = min(max(self.t_tumour + (tumour_influx + tumour_growth - tumour_loss) * dt, 0.0), self.carrying_capacity)

        tumour_volume_l = float(context.get("tumour_volume_l", context.get("V_T", 0.0)))
        tumour_volume_l = max(tumour_volume_l, 1e-9)
        density_raw = min((self.t_tumour * self.tcell_volume_factor) / (tumour_volume_l * 1e6), 1e5)

        if self.smoothing_tau > 0.0:
            alpha = math.exp(-dt / self.smoothing_tau)
            self.density_state = density_raw + (self.density_state - density_raw) * alpha
        else:
            self.density_state = density_raw

        kill_coeff = float(context.get("k_C_T1", 0.0))
        context["tcell_alignment_state"] = self.t_tumour
        context["tcell_density_per_ul"] = max(self.density_state, 0.0)
        context["tcell_kill_rate_hat"] = kill_coeff * self.t_tumour * suppress
