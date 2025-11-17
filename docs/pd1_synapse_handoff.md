# PD‑1 Synapse Handoff

## Background

The QSP surrogate reproduces MATLAB snapshots via the Python runtime (`src/offline`). Tumour observables are already in line, but the PD‑1 checkpoint still depends on the legacy `alignment_driver_block`. We pivoted to a *white-box* driver (`src/offline/modules/pd1_whitebox.py`) so the synapse ODEs and Hill inhibition match the SimBiology model exactly. Details of the broader plan live in `docs/new_alignment_plan.md`, and the semantic-alignment backlog is tracked in `docs/TODO_SEMANTIC_ALIGNMENT.md`.

## Current assets

- MATLAB exporter utilities (under `matlab/scripts/`) can dump arbitrary scenarios via `export_pd1_training_suite.m`. We already generated a 600-scenario training parquet at `artifacts/training/pd1_whitebox_training.parquet`. Each row provides `aPD1`, synapse states, and `H_PD1` trajectories for fitting.
- The Python fitter `scripts/fit_pd1_whitebox.py` reads that parquet, simulates the white-box module, and optimises the synapse parameters. It currently assumes the PD‑1 ODE mirrors MATLAB, but our implementation is still a hybrid approximation.
- Validation scripts:
  - `scripts/dev_pd1_driver_compare.py` plots A-series scenario overlays (A1–A6) for visual inspection.
  - `python -m scripts.validate_surrogate --scenarios A1...` drives the frozen model with `alignment_driver_block`.

## Blocking issues

1. **Synapse ODE parity** – `src/offline/modules/pd1_whitebox.py` still executes a simplified explicit Euler update that diverges from the SimBiology reactions (Reactions 89–92 in `artifacts/matlab_frozen_model/example1/equations.txt`). We need to:
   - Keep state variables in the MATLAB units (molecules per µm²) obtained directly from the snapshot exporter.
   - Integrate the exact RHS (including Reaction 92’s stoichiometry) with a stiff integrator (`solve_ivp(..., method="BDF")`) or the existing segment integrator in `src/offline/simulation.py`.
   - Remove the “scale & clamp” workaround so conservation is enforced by the ODE itself.
2. **Parameter units** – `scripts/fit_pd1_whitebox.py` currently treats `k_PD1_*` and `kd_PD1_*` as if their units were already molar. In `parameters/example1_parameters.json` they are stored as µM or nM. Update the loader so:
   - `k_PD1_*` coefficients are scaled by 1e6 (µM) or 1e9 (nM) before dividing by `d_syn`.
   - Dissociation constants (`kd_PD1_*`) are converted to molarity (µM→1e-6, nM→1e-9), ensuring `koff = kon * KD` matches MATLAB numerically.
   - The cross-arm parameters (`Chi_PD1`, etc.) also honour the nanometer → micrometer conversion.
3. **Regression harness** – once the ODE and units match, re-run:
   - Training fit: `python scripts/fit_pd1_whitebox.py --training-path artifacts/training/pd1_whitebox_training.parquet`.
   - Diagnostics on at least one moderate-dose and one high-dose training scenario by diffing MATLAB vs Python (`artifacts/dev/pd1_compare_*.csv`).
   - A1–A6 validation (`scripts/dev_pd1_driver_compare.py` + `scripts.validate_surrogate`).

## References for onboarding

- `docs/new_alignment_plan.md`: overall alignment strategy and white-box milestones.
- `docs/TODO_SEMANTIC_ALIGNMENT.md`: semantic alignment principles for reproducing SimBiology behaviour.
- `artifacts/matlab_frozen_model/example1/equations.txt`: authoritative definition of the PD‑1 reactions (Reactions 89–92) and `H_PD1` rule.
- `scripts/fit_pd1_whitebox.py`: parameter loader/fitter (needs unit fixes).
- `src/offline/modules/pd1_whitebox.py`: current PD‑1 module (needs full MATLAB ODE + stiff integrator).

## Next steps

你的直觉是对的：**现在应该继续用 `pd1_whitebox.py` 做白盒，但把 M2 里已经踩通坑的 BDF 整套“刚性求解包装”抽出来复用**，而不是在 whitebox 里再造一个小求解器。

下面给你一个可以直接发给新工程师的英文任务说明（包括具体 Step1–4），已经按你写的 Background 风格来写，而且明确要求“复用 M2 的求解器包装”。

---

## PD‑1 Synapse Handoff — Next Engineering Tasks (reuse the M2 stiff solver)

### High‑level goal

Make the PD‑1 white‑box module (`src/offline/modules/pd1_whitebox.py`) behave exactly like the SimBiology synapse ODEs, **without** introducing a completely separate ODE solver. Instead, we want to **reuse the stiff BDF solver configuration that already drives the main frozen model** (the machinery introduced in M2 under `segment_integrator.py` and `simulation.py`). This keeps numerics consistent and makes later refactors (full Python‑side ODEs) much easier.

Concretely:

* The PD‑1 synapse will continue to live as a “white‑box module” in Python.
* The module will call a **shared stiff integration helper** that is factored out of the existing segmented integrator, rather than rolling its own explicit Euler loop.
* Parameter units (kon, KD, etc.) must be cleaned up once, in the same style as the global `units.py` path.

---

### Step 0 — Orientation 

Before touching anything, please skim:

* `docs/new_alignment_plan.md` — overall SimBiology ↔ Python alignment strategy, where PD‑1 is “Step 2”.
* `artifacts/matlab_frozen_model/example1/equations.txt` — authoritative PD‑1 synapse reactions (Reactions 89–92) and the `H_PD1` rule.
* `src/offline/segment_integrator.py` and `src/offline/simulation.py` — how the **main** frozen model uses SciPy’s `solve_ivp` with `SolverConfig`, warm‑start, BDF, max step, etc.
* `src/offline/modules/pd1_whitebox.py` and `scripts/fit_pd1_whitebox.py` — current PD‑1 white‑box implementation and fitter (these still need unit + solver fixes).

You don’t need to understand the whole model; just focus on:

* How the **global** stiff solver is configured (`SolverConfig`, `solve_ivp` usage).
* How the PD‑1 white‑box is currently integrating its 4‑state system.

---

### Step 1 — Extract a reusable stiff ODE helper from the M2 solver

Right now, the stiff integration logic lives inline inside `segment_integrator._warm_start_segment`, `_perform_t0_warm_start`, and `run_segmented_integration`, and is configured via `SolverConfig`.

We want a **small, reusable helper** that PD‑1 (and future white‑box modules) can call, without pulling in events/doses/debouncing.

**1.1. Create a new helper module**

Create a new file, for example:

* `src/offline/stiff_ode.py`

and implement something like:

```python
from typing import Callable, Optional
import numpy as np
from scipy.integrate import solve_ivp

from .segment_integrator import SolverConfig  # reuse existing dataclass

def integrate_local_system(
    rhs: Callable[[float, np.ndarray], np.ndarray],
    y0: np.ndarray,
    t0: float,
    t1: float,
    solver: SolverConfig,
    *,
    max_internal_step_days: float = 1e-4,
) -> np.ndarray:
    """
    Integrate a small stiff ODE system y' = rhs(t, y) from t0 to t1
    using the same BDF configuration as the main frozen-model solver.
    Returns the state at t1.
    """
    span = float(t1 - t0)
    if span <= 0.0:
        return y0.copy()

    # Reuse the same tolerances and method as the main solver.
    max_step = min(float(max_internal_step_days), span)
    sol = solve_ivp(
        rhs,
        (t0, t1),
        y0.copy(),
        method=solver.method,
        rtol=solver.rtol,
        atol=solver.atol,
        max_step=max_step,
        dense_output=False,
        events=None,
        jac_sparsity=None,
        vectorized=False,
    )
    if not sol.success or not sol.y.size:
        # You can refine error handling later
        raise RuntimeError(f"PD1 micro-integration failed: {sol.message}")
    return np.asarray(sol.y[:, -1], dtype=float)
```

Key points:

* Reuse the **same `SolverConfig`** object that `simulate_frozen_model` already builds for the main model.
* Use BDF by default (because `SolverType` is mapped from MATLAB’s `ode15s/ode23s` to `"BDF"`).
* You don’t need any event handling, warm‑start quarantine, or ε‑bump for the PD‑1 micro‑ODE — just short, stiff segments between two macro times.

**1.2. Add unit tests**

Under `tests/`, add e.g. `tests/test_stiff_ode_helper.py`:

* Construct a tiny stiff system (e.g., 2‑state linear system with known analytic solution) and confirm that `integrate_local_system` hits the expected value within tolerance.
* This gives you confidence the helper behaves like the main BDF config without relying on the full frozen model.

---

### Step 2 — Refactor `PD1WhiteboxModel` to call the shared stiff helper

Once `integrate_local_system` exists, the PD‑1 module should **stop** doing its own explicit Euler / sub‑stepping and simply call the shared helper.

**2.1. Implement the exact RHS from equations.txt**

In `src/offline/modules/pd1_whitebox.py`:

* Define a pure function that implements the four‑state RHS exactly as in SimBiology (Reactions 89–92), using **densities** (`molecules/µm²`) and whatever inputs you need (aPD1 concentration, PD‑1/PD‑L totals, etc.).

Example sketch (names to adjust to your current code):

```python
def _pd1_rhs_factory(params, drivers):
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        syn_pd1, syn_pdl1, syn_ab, syn_ab_pd1 = y  # example
        # Compute forward/backward fluxes exactly as in equations.txt
        # Use kon/koff, Chi_PD1, etc. from params, and aPD1(t) from drivers
        # Return np.array([d/dt syn_pd1, d/dt syn_pdl1, ...]) in density units
        ...
    return rhs
```

The goal is that if you plug this RHS into MATLAB’s solver, you would get the same behaviour — no extra scaling or clamping in Python.

**2.2. Make PD1WhiteboxModel a thin wrapper around the helper**

Still in `pd1_whitebox.py`, update `PD1WhiteboxModel.step(...)` (or equivalent) to:

* Build the `rhs` via `_pd1_rhs_factory(...)`.
* Call `integrate_local_system(rhs, y0, t_prev, t_next, solver_config, max_internal_step_days=...)`.
* Store the new state and compute:

  * `H_PD1` (Hill transform of (PD1_PDL1 + PD1_PDL2) densities).
  * The “block fraction” observable used in MATLAB ((PD1_aPD1 + 2·PD1_aPD1_PD1)/total PD1).

Important: **no more hand‑rolled “scale & clamp” logic** — conservation should emerge from the ODE itself.

**2.3. Get a `SolverConfig` into the module**

The PD‑1 module lives below `simulate_frozen_model`, so it doesn’t currently know about `SolverConfig`.

Two simple options:

1. **Pass `solver_config` down via `alignment_driver_block`**

   * `simulate_frozen_model` already creates `solver_config` near the top.
   * When it instantiates `alignment_driver_block`, pass `solver_config` into the module’s constructor or `run(...)` method, and then down into `PD1WhiteboxModel`.

2. **Build a local `SolverConfig` copy with the same values**

   * Less ideal (you risk divergence later), but acceptable if you read the same `SolverOptions` in one place.

Prefer (1) if you can, to keep a single source of truth for solver settings.

---

### Step 3 — Fix PD‑1 parameter units in one place

The solver alone isn’t enough — the **kon/koff/KD units** must match SimBiology. Right now, `scripts/fit_pd1_whitebox.py` and the PD‑1 module treat parameters as if they were already in molar units; the JSON actually stores µM/nM.

The job here is to centralise and clean that up:

**3.1. Extend the parameter loader with explicit unit conversions**

In `scripts/fit_pd1_whitebox.py` (or a dedicated PD‑1 loader), implement something like:

```python
def load_pd1_params(raw_params: Mapping[str, float]) -> PD1Params:
    kon = raw_params["k_PD1_PDL1"] * 1e6  # µM -> M, example
    kd  = raw_params["kd_PD1_PDL1"] * 1e-6  # µM -> M
    koff = kon * kd
    ...
    return PD1Params(kon=kon, koff=koff, ...)
```

* Use the same “unit to M/day” logic you adopted in `units.py` and M2 for other rates.
* Apply nanometer → micrometer scaling for the cross‑arm terms (`Chi_PD1`, etc.) consistently.

**3.2. Pass the cleaned PD‑1 parameters into `PD1WhiteboxModel`**

Don’t let the white‑box module rummage in raw JSON/Snapshot units. Instead:

* `alignment_driver_block` calls `load_pd1_params(...)`.
* It constructs `PD1WhiteboxModel` with a fully unit‑normalised param object.

---

### Step 4 — Regression + validation loop

With (i) the shared stiff helper and (ii) clean PD‑1 units in place, the remaining work is verification; the workflow for the new engineer should be:

**4.1. Single‑scenario diff against MATLAB**

* Use the existing comparison harness (e.g. `pd1_train_0004` in `artifacts/training/pd1_whitebox_training.parquet`).
* For one moderate‑dose and one high‑dose training scenario, generate `artifacts/dev/pd1_compare_*.csv` listing, per time point:

  * `aPD1_molar` (MATLAB)
  * `syn_pd1_*` (all synapse states)
  * `H_PD1` and block fraction (MATLAB vs Python)
  * Pointwise differences

The target is: **curves overlay to within numerical noise**, not just roughly match peaks.

**4.2. Re‑run the PD‑1 training fit**

Once the PD‑1 ODE is pointwise‑aligned for a few scenarios:

```bash
python scripts/fit_pd1_whitebox.py \
  --training-path artifacts/training/pd1_whitebox_training.parquet \
  --output-json artifacts/training/pd1_whitebox_fit.json
```

* This should now refine only small residual differences, not fight a structural mismatch.

**4.3. Validate on A‑series**

Use the standard A‑series overlay tooling:

```bash
python -m scripts.validate_surrogate \
    --scenarios A1 A2 A3 A4 A5 A6 \
    --ic-mode snapshot \
    --module-block alignment_driver_block \
    --dump-flat-debug 5 \
    --param-override alignment_mode=2 --max-rel-err 1e12

python scripts/dev_pd1_driver_compare.py \
    --scenarios A1 A2 A3 A4 A5 A6 \
    --summary-json artifacts/dev/pd1_compare_summary.json
```

* The goal is PD‑1 rel‑RMSE ≪ 0.1 and reasonable peaks across A1–A6.
* Once this is green, you can gradually retire the grey‑box path (`alignment_mode<2`).
