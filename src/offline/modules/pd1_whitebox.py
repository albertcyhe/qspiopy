"""White-box PD-1 checkpoint dynamics helper."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, MutableMapping

import math

SECONDS_PER_DAY = 86400.0


@dataclass
class PD1WhiteboxOutputs:
    """Container for PD-1 white-box step results."""

    occupancy: float
    raw_complexes: float
    complex_density: float
    blocked_fraction: float


@dataclass
class PD1WhiteboxModel:
    total_pd1: float
    total_pdl1: float
    total_pdl2: float
    kon_pd1_pdl1: float
    koff_pd1_pdl1: float
    kon_pd1_pdl2: float
    koff_pd1_pdl2: float
    kon_pd1_ab: float
    koff_pd1_ab: float
    chi_pd1: float
    internalisation: float
    smoothing_tau: float
    occ_smoothed: float
    syn_pd1_pdl1: float
    syn_pd1_pdl2: float
    syn_pd1_ab: float
    syn_pd1_ab_pd1: float
    gamma_c_nivolumab: float
    synapse_area: float
    pd1_50_density: float
    hill_coefficient: float

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

    @classmethod
    def from_context(cls, params: Mapping[str, float], context: Mapping[str, float]) -> "PD1WhiteboxModel":
        def get_context(name: str, fallback: float = 0.0) -> float:
            return float(context.get(name, fallback))

        kon_pd1_pdl1 = float(params.get("kon_PD1_PDL1", 0.0)) * SECONDS_PER_DAY
        kon_pd1_pdl2 = float(params.get("kon_PD1_PDL2", 0.0)) * SECONDS_PER_DAY
        kon_pd1_ab = float(params.get("kon_PD1_aPD1", 0.0)) * SECONDS_PER_DAY
        koff_pd1_pdl1 = float(params.get("koff_PD1_PDL1", 0.0)) * SECONDS_PER_DAY
        koff_pd1_pdl2 = float(params.get("koff_PD1_PDL2", 0.0)) * SECONDS_PER_DAY
        koff_pd1_ab = float(params.get("koff_PD1_aPD1", 0.0)) * SECONDS_PER_DAY
        chi_pd1 = float(params.get("Chi_PD1", 1.0))
        internalisation = float(params.get("pd1_occ_internalization_per_day", 0.0))
        smoothing_tau = float(params.get("pd1_whitebox_tau_days", 0.0))
        gamma_c_nivolumab = float(params.get("gamma_C_nivolumab", 1.0))
        synapse_area = max(float(params.get("A_syn", 1.0)), 1e-6)
        pd1_50_density = max(float(params.get("PD1_50", 1.0)), 1e-12)
        pd1_50_override = params.get("pd1_whitebox_pd1_50_density")
        if pd1_50_override is not None:
            pd1_50_density = max(float(pd1_50_override), 1e-12)
        hill_coefficient = float(params.get("n_PD1", 1.0))

        def _total_density(
            param_key: str,
            context_key: str,
            area_param_key: str,
            default: float = 0.0,
        ) -> float:
            ctx_value = context.get(context_key)
            if ctx_value is not None:
                try:
                    density = float(ctx_value)
                except (TypeError, ValueError):
                    density = float("nan")
                if math.isfinite(density) and density > 0.0:
                    return density
            param_value = params.get(param_key, default)
            if param_value is None:
                return max(default, 0.0)
            try:
                molecules = float(param_value)
            except (TypeError, ValueError):
                molecules = default
            area_value = params.get(area_param_key, synapse_area)
            try:
                area = float(area_value)
            except (TypeError, ValueError):
                area = synapse_area
            if area <= 0.0:
                area = synapse_area
            return max(molecules / area, 0.0)

        total_pd1 = _total_density("T_PD1_total", "syn_pd1_total", "A_Tcell")
        total_pdl1 = _total_density("C_PDL1_total", "syn_pdl1_total", "A_cell")
        total_pdl2 = _total_density("C_PDL2_total", "syn_pdl2_total", "A_cell")

        syn_pd1_pdl1 = get_context("syn_pd1_pdl1", 0.0)
        syn_pd1_pdl2 = get_context("syn_pd1_pdl2", 0.0)
        syn_pd1_ab = get_context("syn_pd1_apd1", 0.0)
        syn_pd1_ab_pd1 = get_context("syn_pd1_apd1_pd1", 0.0)
        initial_raw = synapse_area * max(syn_pd1_pdl1, 0.0)
        occ_smoothed = cls._hill_value(initial_raw, pd1_50_density, hill_coefficient)

        return cls(
            total_pd1=max(total_pd1, 0.0),
            total_pdl1=max(total_pdl1, 0.0),
            total_pdl2=max(total_pdl2, 0.0),
            kon_pd1_pdl1=kon_pd1_pdl1,
            koff_pd1_pdl1=koff_pd1_pdl1,
            kon_pd1_pdl2=kon_pd1_pdl2,
            koff_pd1_pdl2=koff_pd1_pdl2,
            kon_pd1_ab=kon_pd1_ab,
            koff_pd1_ab=koff_pd1_ab,
            chi_pd1=chi_pd1,
            internalisation=max(internalisation, 0.0),
            smoothing_tau=max(smoothing_tau, 0.0),
            occ_smoothed=occ_smoothed,
            syn_pd1_pdl1=max(syn_pd1_pdl1, 0.0),
            syn_pd1_pdl2=max(syn_pd1_pdl2, 0.0),
            syn_pd1_ab=max(syn_pd1_ab, 0.0),
            syn_pd1_ab_pd1=max(syn_pd1_ab_pd1, 0.0),
            synapse_area=synapse_area,
            pd1_50_density=pd1_50_density,
            hill_coefficient=hill_coefficient,
            gamma_c_nivolumab=max(gamma_c_nivolumab, 1e-6),
        )

    def _pd1_free(self) -> float:
        bound = self.syn_pd1_pdl1 + self.syn_pd1_pdl2 + self.syn_pd1_ab + 2.0 * self.syn_pd1_ab_pd1
        return max(self.total_pd1 - bound, 0.0)

    def _pdl1_free(self) -> float:
        return max(self.total_pdl1 - self.syn_pd1_pdl1, 0.0)

    def _pdl2_free(self) -> float:
        return max(self.total_pdl2 - self.syn_pd1_pdl2, 0.0)

    def step(self, antibody_concentration_molar: float, delta_t: float) -> PD1WhiteboxOutputs:
        delta_t = max(float(delta_t), 0.0)
        if delta_t <= 0.0:
            raw_density = max(self.syn_pd1_pdl1 + self.syn_pd1_pdl2, 0.0)
            blocked = self._blocked_fraction()
            return PD1WhiteboxOutputs(self.occ_smoothed, raw_density, self.syn_pd1_pdl1, blocked_fraction=blocked)

        pd1_free = self._pd1_free()
        pdl1_free = self._pdl1_free()
        pdl2_free = self._pdl2_free()
        ab_concentration = max(float(antibody_concentration_molar), 0.0)
        gamma_c = self.gamma_c_nivolumab
        ab_effective = ab_concentration / gamma_c

        d_pd1_pdl1 = (
            self.kon_pd1_pdl1 * pd1_free * pdl1_free
            - self.koff_pd1_pdl1 * self.syn_pd1_pdl1
            - self.internalisation * self.syn_pd1_pdl1
        )
        d_pd1_pdl2 = (
            self.kon_pd1_pdl2 * pd1_free * pdl2_free
            - self.koff_pd1_pdl2 * self.syn_pd1_pdl2
            - self.internalisation * self.syn_pd1_pdl2
        )
        d_pd1_ab = (
            2.0 * self.kon_pd1_ab * pd1_free * ab_effective
            - self.koff_pd1_ab * self.syn_pd1_ab
            - self.internalisation * self.syn_pd1_ab
        )
        d_pd1_ab_pd1 = (
            self.chi_pd1 * self.kon_pd1_ab * pd1_free * self.syn_pd1_ab
            - 2.0 * self.koff_pd1_ab * self.syn_pd1_ab_pd1
            - self.internalisation * self.syn_pd1_ab_pd1
        )

        self.syn_pd1_pdl1 = min(max(self.syn_pd1_pdl1 + d_pd1_pdl1 * delta_t, 0.0), self.total_pd1)
        self.syn_pd1_pdl2 = min(max(self.syn_pd1_pdl2 + d_pd1_pdl2 * delta_t, 0.0), self.total_pd1)
        self.syn_pd1_ab = min(max(self.syn_pd1_ab + d_pd1_ab * delta_t, 0.0), self.total_pd1)
        self.syn_pd1_ab_pd1 = min(
            max(self.syn_pd1_ab_pd1 + d_pd1_ab_pd1 * delta_t, 0.0),
            self.total_pd1,
        )
        bound_pd1 = (
            self.syn_pd1_pdl1
            + self.syn_pd1_pdl2
            + self.syn_pd1_ab
            + 2.0 * self.syn_pd1_ab_pd1
        )
        if bound_pd1 > self.total_pd1 and bound_pd1 > 0.0:
            scale = self.total_pd1 / bound_pd1
            self.syn_pd1_pdl1 *= scale
            self.syn_pd1_pdl2 *= scale
            self.syn_pd1_ab *= scale
            self.syn_pd1_ab_pd1 *= scale

        raw_density = max(self.syn_pd1_pdl1 + self.syn_pd1_pdl2, 0.0)
        raw_occ = self._hill_value(raw_density, self.pd1_50_density, self.hill_coefficient)
        raw_occ = min(max(raw_occ, 0.0), 1.0)
        if self.smoothing_tau > 0.0:
            alpha = math.exp(-delta_t / self.smoothing_tau)
            self.occ_smoothed = raw_occ + (self.occ_smoothed - raw_occ) * alpha
        else:
            self.occ_smoothed = raw_occ
        blocked = self._blocked_fraction()
        return PD1WhiteboxOutputs(self.occ_smoothed, raw_density, self.syn_pd1_pdl1, blocked_fraction=blocked)

    def _blocked_fraction(self) -> float:
        if self.total_pd1 <= 0.0:
            return 0.0
        blocked = (self.syn_pd1_ab + 2.0 * self.syn_pd1_ab_pd1) / self.total_pd1
        return min(max(blocked, 0.0), 1.0)

    def writeback(self, context: MutableMapping[str, float]) -> None:
        context["syn_pd1_total"] = self.total_pd1
        context["syn_pdl1_total"] = self.total_pdl1
        context["syn_pdl2_total"] = self.total_pdl2
        context["syn_pd1_pdl1"] = self.syn_pd1_pdl1
        context["syn_pd1_pdl2"] = self.syn_pd1_pdl2
        context["syn_pd1_apd1"] = self.syn_pd1_ab
        context["syn_pd1_apd1_pd1"] = self.syn_pd1_ab_pd1
