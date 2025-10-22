"""Simplified reproductions of the QSP-IO tutorial case studies.

The original MATLAB tutorial assembles a rich SimBiology model to demonstrate
checkpoint inhibition scenarios.  Re-implementing every mechanistic detail is
well beyond the scope of the offline Python surrogate; however, having a
lightweight approximation is invaluable for regression testing and for
exercising the parameter catalogues that ship with the publication (as well as
external catalogues such as the TNBC extension project).

The routines below operate on top of :mod:`src.offline.parameter_loader` and
simulate a reduced three-state system consisting of viable cancer cells, dead
cancer debris, and activated tumour-infiltrating T cells.  While compact, the
model honours several qualitative behaviours highlighted in the tutorial:

* Logistic cancer growth capped by the carrying capacity described in the
  parameter catalogue.
* Cytotoxic suppression proportional to the availability of effector T cells
  and modulated by PD-1/PD-L1 checkpoint inhibition.
* Dead-cell clearance contributing to the pseudo-progression phenomenon when
  calculating tumour volume.

The approximation is certainly cruder than the full SimBiology pipeline, but
it provides a deterministic and reproducible scaffold for experimenting with
the tutorial parameterisations directly from Python.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Dict, Iterable, Literal

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from . import parameter_loader


@dataclass(frozen=True)
class ScenarioResult:
    """Container describing the output of a simulation run."""

    time_days: np.ndarray
    cancer_cells: np.ndarray
    dead_cells: np.ndarray
    t_cells: np.ndarray
    tumour_volume_l: np.ndarray
    tumour_diameter_cm: np.ndarray
    pd1_occupancy: np.ndarray
    tcell_density_per_ul: np.ndarray

    def to_frame(self) -> pd.DataFrame:
        """Return the simulation as a tidy :class:`pandas.DataFrame`."""

        return pd.DataFrame(
            {
                "time_days": self.time_days,
                "cancer_cells": self.cancer_cells,
                "dead_cells": self.dead_cells,
                "t_cells": self.t_cells,
                "tumour_volume_l": self.tumour_volume_l,
                "tumour_diameter_cm": self.tumour_diameter_cm,
                "pd1_occupancy": self.pd1_occupancy,
                "tcell_density_per_ul": self.tcell_density_per_ul,
            }
        )


def _checkpoint_inhibition(params: Dict[str, float], *, therapy: bool) -> float:
    """Approximate the PD-1 mediated inhibition factor.

    The SimBiology model explicitly simulates the PD-1 / PD-L1 binding kinetics
    inside the immune synapse.  To keep the Python surrogate tractable we
    collapse those dynamics into a single Hill-type term which depends on the
    number of complexes (approximated by the geometric mean of the available
    ligands) and optionally reduced if anti PD-1 therapy is active.
    """

    pd1_50 = params.get("PD1_50", 1.0)
    n_pd1 = params.get("n_PD1", 1.0)
    pd1_total = params.get("T_PD1_total", 1e4)
    pdl1_total = params.get("C_PDL1_total", 1e4)

    complexes = math.sqrt(pd1_total * pdl1_total)
    inhibition = complexes ** n_pd1 / (complexes ** n_pd1 + pd1_50 ** n_pd1)

    if therapy:
        kd = params.get("kd_PD1_aPD1", 1.0)
        chi = params.get("Chi_PD1", 1.0)
        inhibition *= 1.0 / (1.0 + chi / max(kd, 1e-9))
    return float(max(min(inhibition, 1.0), 0.0))


def simulate_tutorial(
    parameter_paths: Iterable[Path | str],
    *,
    days: float = 400.0,
    therapy: Literal["none", "anti_pd1"] = "anti_pd1",
) -> ScenarioResult:
    """Simulate the reduced-order tutorial model.

    Parameters
    ----------
    parameter_paths:
        One or more JSON catalogues.  Later entries override earlier ones –
        this mirrors the MATLAB behaviour where project-specific files extend
        the base catalogue.
    days:
        Length of the simulation horizon.
    therapy:
        ``"none"`` keeps the checkpoint inhibition factor untouched whereas
        ``"anti_pd1"`` applies the nivolumab-driven reduction discussed in the
        tutorial examples.
    """

    params = parameter_loader.combine_parameter_sets(parameter_paths)

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

    inhibition = _checkpoint_inhibition(params, therapy=therapy == "anti_pd1")

    def rhs(_t: float, state: np.ndarray) -> np.ndarray:
        cancer, dead, t_cells = state
        cytotoxic_term = k_tcell * t_cells / (t_cells + cancer + 1e-9) * (1.0 - inhibition)

        d_cancer = (
            k_growth * cancer * (1.0 - cancer / max(c_max, 1e-9))
            - (k_innate + cytotoxic_term) * cancer
        )
        d_dead = (k_innate + cytotoxic_term) * cancer - k_clear * dead
        recruitment = params.get("Q_nT_thym", 0.0)
        carrying = params.get("div", 1.0)
        d_t_cells = (
            recruitment
            + k_pro * t_cells * (1.0 - t_cells / max(carrying, 1.0))
            - k_death * t_cells
        )
        return np.array([d_cancer, d_dead, d_t_cells], dtype=float)

    time_span = (0.0, float(days))
    time_eval = np.linspace(time_span[0], time_span[1], 401)

    sol = solve_ivp(
        rhs,
        time_span,
        np.array([initial_cancer, 0.0, t_cells0], dtype=float),
        t_eval=time_eval,
        method="LSODA",
        vectorized=False,
        rtol=1e-6,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(f"Integration failed: {sol.message}")

    cancer = np.maximum(sol.y[0], 0.0)
    dead = np.maximum(sol.y[1], 0.0)
    t_cells = np.maximum(sol.y[2], 0.0)

    volumes = parameter_loader.tumour_volume_liters(
        cancer, params, t_cells=t_cells, dead_cells=dead
    )
    diameters = parameter_loader.tumour_diameter_from_volume(volumes)
    pd1_series = np.full_like(sol.t, inhibition, dtype=float)
    density_per_ul = t_cells / np.maximum(volumes * 1e6, 1e-12)

    return ScenarioResult(
        time_days=sol.t,
        cancer_cells=cancer,
        dead_cells=dead,
        t_cells=t_cells,
        tumour_volume_l=volumes,
        tumour_diameter_cm=diameters,
        pd1_occupancy=pd1_series,
        tcell_density_per_ul=density_per_ul,
    )

