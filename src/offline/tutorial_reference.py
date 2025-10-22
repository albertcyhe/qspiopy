"""Reference implementation mirroring the MATLAB tutorial equations.

The purpose of this module is to provide a high-fidelity baseline that follows
the dynamics outlined in the original QSP-IO tutorial article.  While it still
operates entirely in Python, the structure mirrors the published differential
equations and uses a relatively stiff solver setup to mimic the behaviour of
the SimBiology workflow.  The resulting trajectories act as a quantitative
reference when validating the lighter-weight surrogate implemented in
``tutorial_sim``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import numpy as np
from scipy.integrate import solve_ivp

from . import parameter_loader
from .tutorial_sim import ScenarioResult, _checkpoint_inhibition


@dataclass(frozen=True)
class ReferenceSettings:
    """Solver configuration approximating the SimBiology defaults."""

    rtol: float = 1e-7
    atol: float = 1e-10
    method: str = "BDF"
    dense_output: bool = True


def simulate_reference(
    parameter_paths: Iterable[Path | str],
    *,
    days: float = 400.0,
    therapy: Literal["none", "anti_pd1"] = "anti_pd1",
    settings: ReferenceSettings | None = None,
) -> ScenarioResult:
    """Simulate the mechanistic reference model for the tutorial scenarios."""

    params = parameter_loader.combine_parameter_sets(parameter_paths)
    settings = settings or ReferenceSettings()

    k_growth = params["k_C_growth"]
    c_max = params["C_max"]
    k_innate = params.get("k_C_death", 0.0)
    k_tcell = params.get("k_C_Tcell", 0.0)
    k_clear = params.get("k_cell_clear", 0.01)
    k_pro = params.get("k_pro", 0.5)
    k_death = params.get("k_death", 0.5)

    initial_cancer = parameter_loader.estimate_initial_tumour_cells(params)
    initial_volume = parameter_loader.tumour_volume_liters(initial_cancer, params)
    rho_nT = params.get("rho_nT", 5.0e5)
    t_cells0 = rho_nT * (initial_volume * 1e3)

    base_inhibition = _checkpoint_inhibition(params, therapy=False)
    therapy_inhibition = _checkpoint_inhibition(params, therapy=True)

    tau_pd1 = params.get("tau_PD1", 5.0)
    activation_gain = params.get("activation_gain", 0.0)

    target_inhibition = therapy_inhibition if therapy == "anti_pd1" else base_inhibition

    def rhs(_t: float, state: np.ndarray) -> np.ndarray:
        cancer, dead, t_cells, occ = state

        effective_occ = max(min(occ, 1.0), 0.0)
        cytotoxic_term = k_tcell * t_cells / (t_cells + cancer + 1e-9) * (1.0 - effective_occ)

        d_cancer = (
            k_growth * cancer * (1.0 - cancer / max(c_max, 1e-9))
            - (k_innate + cytotoxic_term) * cancer
        )
        d_dead = (k_innate + cytotoxic_term) * cancer - k_clear * dead

        recruitment = params.get("Q_nT_thym", 0.0)
        carrying = params.get("div", 1.0)
        activation_boost = activation_gain * (1.0 - effective_occ)
        d_t_cells = (
            recruitment
            + k_pro * t_cells * (1.0 - t_cells / max(carrying, 1.0))
            + activation_boost * cancer
            - k_death * t_cells
        )
        d_occ = (target_inhibition - occ) / max(tau_pd1, 1e-6)

        return np.array([d_cancer, d_dead, d_t_cells, d_occ], dtype=float)

    time_span = (0.0, float(days))
    time_eval = np.linspace(time_span[0], time_span[1], 401)

    sol = solve_ivp(
        rhs,
        time_span,
        np.array([initial_cancer, 0.0, t_cells0, base_inhibition], dtype=float),
        t_eval=time_eval,
        method=settings.method,
        rtol=settings.rtol,
        atol=settings.atol,
        dense_output=settings.dense_output,
    )
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")

    cancer = np.maximum(sol.y[0], 0.0)
    dead = np.maximum(sol.y[1], 0.0)
    t_cells = np.maximum(sol.y[2], 0.0)
    occupancy = np.clip(sol.y[3], 0.0, 1.0)

    volumes = parameter_loader.tumour_volume_liters(
        cancer, params, t_cells=t_cells, dead_cells=dead
    )
    diameters = parameter_loader.tumour_diameter_from_volume(volumes)
    density_per_ul = t_cells / np.maximum(volumes * 1e6, 1e-12)

    return ScenarioResult(
        time_days=sol.t,
        cancer_cells=cancer,
        dead_cells=dead,
        t_cells=t_cells,
        tumour_volume_l=volumes,
        tumour_diameter_cm=diameters,
        pd1_occupancy=occupancy,
        tcell_density_per_ul=density_per_ul,
    )

