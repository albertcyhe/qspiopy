"""White-box PD-1 checkpoint dynamics helper."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, MutableMapping, Optional

import logging
import math
import numpy as np

from ..errors import NumericsError
from ..segment_integrator import SolverConfig
from ..stiff_ode import integrate_local_system
from .pd1_params import PD1Params

logger = logging.getLogger(__name__)


def pd1_synapse_reaction_terms(
    y: np.ndarray,
    params: PD1Params,
    *,
    total_pd1_density: float,
    total_pdl1_density: float,
    total_pdl2_density: float,
    ab_effective_molar: float,
) -> tuple[float, float, float, float, float, float, float]:
    """Return individual reaction fluxes and free species densities."""
    pdl1, pdl2, ab, ab_pd1 = y
    bound = pdl1 + pdl2 + ab + 2.0 * ab_pd1
    pd1_free = max(total_pd1_density - bound, 0.0)
    pdl1_free = max(total_pdl1_density - pdl1, 0.0)
    pdl2_free = max(total_pdl2_density - pdl2, 0.0)

    reaction_89 = params.kon_pd1_pdl1 * pd1_free * pdl1_free - params.koff_pd1_pdl1 * pdl1
    reaction_90 = params.kon_pd1_pdl2 * pd1_free * pdl2_free - params.koff_pd1_pdl2 * pdl2
    reaction_91 = 2.0 * params.kon_pd1_ab * pd1_free * ab_effective_molar - params.koff_pd1_ab * ab
    reaction_92 = params.chi_pd1 * params.kon_pd1_ab * pd1_free * ab - 2.0 * params.koff_pd1_ab * ab_pd1
    return reaction_89, reaction_90, reaction_91, reaction_92, pd1_free, pdl1_free, pdl2_free


def pd1_synapse_rhs(
    y: np.ndarray,
    params: PD1Params,
    *,
    total_pd1_density: float,
    total_pdl1_density: float,
    total_pdl2_density: float,
    ab_effective_molar: float,
) -> np.ndarray:
    """Evaluate RHS for the PD-1 synapse system in density units."""
    reaction_89, reaction_90, reaction_91, reaction_92, _, _, _ = pd1_synapse_reaction_terms(
        y,
        params,
        total_pd1_density=total_pd1_density,
        total_pdl1_density=total_pdl1_density,
        total_pdl2_density=total_pdl2_density,
        ab_effective_molar=ab_effective_molar,
    )
    pdl1, pdl2, ab, ab_pd1 = y
    dpdl1 = reaction_89 - params.internalisation_per_day * pdl1
    dpdl2 = reaction_90 - params.internalisation_per_day * pdl2
    dab = reaction_91 - reaction_92 - params.internalisation_per_day * ab
    dab_pd1 = reaction_92 - params.internalisation_per_day * ab_pd1
    return np.array([dpdl1, dpdl2, dab, dab_pd1], dtype=float)


@dataclass
class PD1WhiteboxOutputs:
    occupancy: float
    raw_complexes: float
    complex_density: float
    blocked_fraction: float


@dataclass
class PD1WhiteboxModel:
    params: PD1Params
    solver_config: SolverConfig
    total_pd1_density: float
    total_pdl1_density: float
    total_pdl2_density: float
    occ_smoothed: float
    syn_pd1_pdl1: float
    syn_pd1_pdl2: float
    syn_pd1_ab: float
    syn_pd1_ab_pd1: float
    time_days: float
    _debug_enabled: bool = False

    @classmethod
    def from_context(
        cls,
        params: PD1Params,
        context: Mapping[str, float],
        solver_config: SolverConfig | None = None,
    ) -> PD1WhiteboxModel:
        def _density(key: str, fallback: float) -> float:
            value = context.get(key)
            if value is None:
                return fallback
            try:
                return max(float(value), 0.0)
            except (TypeError, ValueError):
                return fallback

        pd1_density = _density("syn_pd1_total", params.pd1_surface_density())
        pdl1_density = _density("syn_pdl1_total", params.pdl1_surface_density())
        pdl2_density = _density("syn_pdl2_total", params.pdl2_surface_density())
        syn_pd1_pdl1 = _density("syn_pd1_pdl1", 0.0)
        syn_pd1_pdl2 = _density("syn_pd1_pdl2", 0.0)
        syn_pd1_ab = _density("syn_pd1_apd1", 0.0)
        syn_pd1_ab_pd1 = _density("syn_pd1_apd1_pd1", 0.0)
        initial_raw = max(syn_pd1_pdl1 + syn_pd1_pdl2, 0.0)
        occ_smoothed = cls._hill_value(initial_raw, params.pd1_50_density, params.hill_coefficient)
        time_days = float(context.get("time_days", context.get("t", 0.0)))
        cfg = solver_config or SolverConfig(
            method="BDF",
            rtol=params.solver_rtol,
            atol=params.solver_atol,
            max_step=params.max_step_days,
            seed=None,
        )
        return cls(
            params=params,
            solver_config=cfg,
            total_pd1_density=pd1_density,
            total_pdl1_density=pdl1_density,
            total_pdl2_density=pdl2_density,
            occ_smoothed=occ_smoothed,
            syn_pd1_pdl1=syn_pd1_pdl1,
            syn_pd1_pdl2=syn_pd1_pdl2,
            syn_pd1_ab=syn_pd1_ab,
            syn_pd1_ab_pd1=syn_pd1_ab_pd1,
            time_days=time_days,
        )

    @staticmethod
    def _hill_value(raw_complexes: float, half_max: float, hill_coefficient: float) -> float:
        signal = max(float(raw_complexes), 0.0)
        half = max(float(half_max), 1e-12)
        n = max(float(hill_coefficient), 1e-6)
        numerator = math.pow(signal, n)
        denominator = numerator + math.pow(half, n)
        if denominator <= 0.0:
            return 0.0
        return numerator / denominator

    def _pd1_free(self) -> float:
        bound = self.syn_pd1_pdl1 + self.syn_pd1_pdl2 + self.syn_pd1_ab + 2.0 * self.syn_pd1_ab_pd1
        return max(self.total_pd1_density - bound, 0.0)

    def _pdl1_free(self) -> float:
        return max(self.total_pdl1_density - self.syn_pd1_pdl1, 0.0)

    def _pdl2_free(self) -> float:
        return max(self.total_pdl2_density - self.syn_pd1_pdl2, 0.0)

    def step(self, antibody_concentration_molar: float, delta_t: float) -> PD1WhiteboxOutputs:
        delta_t = max(float(delta_t), 0.0)
        if delta_t <= 0.0:
            blocked = self._blocked_fraction()
            raw_density = max(self.syn_pd1_pdl1 + self.syn_pd1_pdl2, 0.0)
            return PD1WhiteboxOutputs(self.occ_smoothed, raw_density, self.syn_pd1_pdl1, blocked_fraction=blocked)

        start_time = self.time_days
        end_time = start_time + delta_t
        self._integrate_segment(antibody_concentration_molar, start_time, end_time)
        self.time_days = end_time

        raw_density = max(self.syn_pd1_pdl1 + self.syn_pd1_pdl2, 0.0)
        raw_occ = self._hill_value(raw_density, self.params.pd1_50_density, self.params.hill_coefficient)
        raw_occ = min(max(raw_occ, 0.0), 1.0)
        if self.params.smoothing_tau_days > 0.0:
            alpha = math.exp(-delta_t / self.params.smoothing_tau_days)
            self.occ_smoothed = raw_occ + (self.occ_smoothed - raw_occ) * alpha
        else:
            self.occ_smoothed = raw_occ
        blocked = self._blocked_fraction()
        return PD1WhiteboxOutputs(self.occ_smoothed, raw_density, self.syn_pd1_pdl1, blocked_fraction=blocked)

    def _integrate_segment(self, antibody_concentration_molar: float, start_time: float, end_time: float) -> None:
        if end_time <= start_time:
            return

        ab_concentration = max(float(antibody_concentration_molar), 0.0)
        ab_effective = ab_concentration / max(self.params.gamma_c_nivolumab, 1e-9)

        y0 = np.array(
            [self.syn_pd1_pdl1, self.syn_pd1_pdl2, self.syn_pd1_ab, self.syn_pd1_ab_pd1],
            dtype=float,
        )

        def rhs(_, y: np.ndarray) -> np.ndarray:
            return pd1_synapse_rhs(
                y,
                self.params,
                total_pd1_density=self.total_pd1_density,
                total_pdl1_density=self.total_pdl1_density,
                total_pdl2_density=self.total_pdl2_density,
                ab_effective_molar=ab_effective,
            )

        debug_cb = self._solver_debug if self._debug_enabled else None

        def _advance(t0: float, t1: float, state: np.ndarray, depth: int = 0) -> np.ndarray:
            span = t1 - t0
            if span <= 1e-12:
                return state
            try:
                internal_step = self.params.max_step_days
                if span > 0.0:
                    internal_step = min(internal_step, max(span * 0.25, 1e-9))
                return integrate_local_system(
                    rhs,
                    state,
                    t0,
                    t1,
                    solver=self.solver_config,
                    max_internal_step_days=internal_step,
                    debug=debug_cb,
                )
            except NumericsError:
                if depth >= 6 or span <= 1e-6:
                    raise
                mid = t0 + 0.5 * span
                mid_state = _advance(t0, mid, state, depth + 1)
                return _advance(mid, t1, mid_state, depth + 1)

        final = _advance(start_time, end_time, y0)
        final = np.maximum(final, 0.0)
        self.syn_pd1_pdl1 = float(final[0])
        self.syn_pd1_pdl2 = float(final[1])
        self.syn_pd1_ab = float(final[2])
        self.syn_pd1_ab_pd1 = float(final[3])

    def _blocked_fraction(self) -> float:
        if self.total_pd1_density <= 0.0:
            return 0.0
        blocked = (self.syn_pd1_ab + 2.0 * self.syn_pd1_ab_pd1) / self.total_pd1_density
        return min(max(blocked, 0.0), 1.0)

    def writeback(self, context: MutableMapping[str, float]) -> None:
        # Update both canonical aliases and the original snapshot identifiers so that
        # FrozenModel.sync_state_from_context persists the modified synapse states.
        context["syn_T1_C1.PD1"] = self.total_pd1_density
        context["syn_pd1_total"] = self.total_pd1_density
        context["syn_T1_C1.PDL1"] = self.total_pdl1_density
        context["syn_pdl1_total"] = self.total_pdl1_density
        context["syn_T1_C1.PDL2"] = self.total_pdl2_density
        context["syn_pdl2_total"] = self.total_pdl2_density
        context["syn_T1_C1.PD1_PDL1"] = self.syn_pd1_pdl1
        context["syn_pd1_pdl1"] = self.syn_pd1_pdl1
        context["syn_T1_C1.PD1_PDL2"] = self.syn_pd1_pdl2
        context["syn_pd1_pdl2"] = self.syn_pd1_pdl2
        context["syn_T1_C1.PD1_aPD1"] = self.syn_pd1_ab
        context["syn_pd1_apd1"] = self.syn_pd1_ab
        context["syn_T1_C1.PD1_aPD1_PD1"] = self.syn_pd1_ab_pd1
        context["syn_pd1_apd1_pd1"] = self.syn_pd1_ab_pd1

    def _solver_debug(self, payload: Mapping[str, object]) -> None:
        if not self._debug_enabled:
            return
        logger.debug(
            "pd1_whitebox solver span=(%.6g, %.6g) status=%s message=%s nfev=%s njev=%s nlu=%s",
            payload.get("t0"),
            payload.get("t1"),
            payload.get("status"),
            payload.get("message"),
            payload.get("nfev"),
            payload.get("njev"),
            payload.get("nlu"),
        )
