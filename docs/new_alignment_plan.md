# Snapshot Alignment Plan (Updated 2025‑11‑13)

This file tracks the high-level plan for bringing the MATLAB snapshots and the Python surrogate back into numeric agreement. For a full handoff narrative see `docs/project_handoff.md`.

---

## Milestone Snapshot

| Milestone | Status | Notes |
| --- | --- | --- |
| **M1 – Event semantics** | ✅ | Segmented integration, trigger specs, pre/post logging. |
| **M2 – Units & parameters** | ✅ | `units.py` canonicalises to day/L/mol/M; dose application goes through `normalise_dose_to_species`. |
| **M3 – Snapshot runtime** | ✅ | Snapshot/target ICs, warm-start quarantine, module hooks, `--dump-flat-debug`. |
| **M4 – Dynamic volume & PD‑1 alignment** | 🚧 | `alignment_driver_block` implements surrogate PK/PD‑1/geometry but needs white-box parity + calibration. |
| **M5 – Numeric gates & CI** | 🚧 | `validate_surrogate` stable; gates still red for A-series. |

---

## Current State (A-series focus)

* Snapshot pipeline is solid: `FrozenModel` loads full CSV/TXT artifacts, the integrator reproduces example1/example2 baselines, and instrumentation shows all intermediate signals.
* PD‑1 / geometry outputs are still driven by runtime modules:
  - `alignment_driver_block` now runs a dose-driven PK state, PD‑1 binding ODE (kon/koff/k_int), and a logistic tumour-volume follower; parameters live in both `parameters/example1_parameters.json` and the exported snapshot CSV.
  - Default bridge/geometry modules are disabled for scenario A1 to keep PD‑1/tumour volume under a single entry point.
* MATLAB side `show_alignment_drivers.m` vs Python `scripts/dump_snapshot_semantics.py` outputs are stored under `artifacts/show_alignment_example1.txt` and `artifacts/show_snapshot_example1_filtered.txt`. They agree on the equations we care about, but exporter coverage still needs to be double‑checked before turning off the alignment layer.

---

## Open Issues

1. **Exporter parity** — verify the MATLAB exporter writes every reaction/rule/event to `reactions.csv`/`equations.txt`. If not, fix the exporter and re-freeze the snapshots.
2. **Runtime overrides** — ensure A-series scenarios only activate `alignment_driver_block` so white-box rules (e.g. `H_PD1_C1`, `V_T`) aren’t disabled.
3. **Parameter calibration** — alignment driver parameters (`pd1_occ_*`, `pd1_pk_*`, `geom_*`) are still defaults; rel‑L2 for `pd1_occupancy`, `tumour_volume_l`, and `tcell_density_per_ul` remains O(1).

---

## Next Actions

1. **Exporter audit**  
   - Run `matlab/scripts/show_alignment_drivers.m` and check the corresponding snapshot dump via  
     `scripts/dump_snapshot_semantics.py artifacts/matlab_frozen_model/example1 --keyword PD1 --keyword V_T --keyword nivol`.  
   - If any reaction/rule is missing, update the MATLAB exporter and re-export `artifacts/matlab_frozen_model/example1`.

2. **Alignment driver calibration**  
   - Keep CLI usage to `python -m scripts.validate_surrogate --scenarios A1 --ic-mode snapshot --module-block alignment_driver_block --dump-flat-debug 20 --numeric-gates`.  
   - Use `scripts.fit_observables.py` (two-stage fitting: PD‑1 first, geometry second) to tune the new parameters and write them back to `parameters/*.json` + snapshot CSV.

3. **Scale to other scenarios**  
   - Once A1 gates pass, run the same driver on A2–A6/B. Introduce scenario-specific overrides only when the MATLAB references truly diverge.

4. **Retire grey-box modules**  
   - When exporter parity is confirmed and the alignment driver reproduces MATLAB shapes, mark `alignment_driver_block` as “debug only” and revert to pure snapshot semantics.

---

## Useful Commands

```bash
# Dump snapshot semantics for PD‑1 / geometry keywords
scripts/dump_snapshot_semantics.py artifacts/matlab_frozen_model/example1 \
  --keyword PD1 --keyword nivol --keyword V_T

# Run surrogate with alignment driver and diagnostics
python -m scripts.validate_surrogate --scenarios A1 --ic-mode snapshot \
  --module-block alignment_driver_block --dump-flat-debug 20 --numeric-gates

# Fit alignment-driver parameters
python -m scripts.fit_observables --scenario A1 \
  --module-block alignment_driver_block --param pd1_occ_kon_scale=1e-8:1e-4 \
  --param geom_growth_per_day=0.005:0.05 ...
```

Keep this plan lean—detailed status, risks, and historical notes live in `docs/project_handoff.md` and `docs/m4_geometry_status.md`.

---

## Immediate Execution Plan (PD‑1 & T cell focus)

### Step 0 – Freeze the baseline *(completed)*

```
python -m scripts.validate_surrogate \
  --scenarios A1 A2 A3 A4 A5 A6 \
  --ic-mode snapshot \
  --module-block alignment_driver_block \
  --dump-flat-debug 5
```

Acceptance criteria:
1. All six A-series scenarios complete deterministically with the same command-line flags.
2. `artifacts/validation/metrics.csv` shows PK/plasma observables and any non‑PD‑1/T-cell outputs comfortably below rel_RMSE ~= 1 (ideally ≪1).

Outcome: the “unknowns” shrink to PD‑1 occupancy and T-cell geometry only.

### Step 1 – Multi-scenario PD‑1 driver comparison *(completed)*

Purpose: understand whether the alignment driver itself (not just `pd1_occupancy`) is saturating or has the wrong time scale.

Actions:
- Re-run Step 0 with `--dump-flat-debug 50` so `A*_flat_debug_*.csv` captures `pd1_alignment_*` signals.
- Write a small plotting script/notebook that, for each A1–A6 scenario, overlays:
  * Surrogate PK state (`pd1_alignment_pk_state`)  
  * Surrogate “effective binding signal” (e.g. `pd1_filter_input` / `pd1_alignment_concentration_M`)  
  * Surrogate `pd1_occupancy`  
  * MATLAB reference `pd1_occupancy`
- Answer:
  * **A:** Are high-dose scenarios (A3/A5/A6) pinned at the driver ceiling immediately after the first dose?
  * **B:** Is the surrogate settling on a drastically different time constant vs MATLAB?
  * **C:** Do different A-regimens yield distinct driver waveforms, or is the grey-box structure incapable of differentiating dose ladders?

Result informs whether PD‑1 needs more poles/turnover mechanics before deeper white-box work.

### Step 2 – Partial white-boxing of the PD‑1 block *(scaffolded; pausing)*

1. Keep exporting PD‑1 reactions/rules with `show_alignment_drivers.m` and `scripts/dump_snapshot_semantics.py`.
2. Enable `pd1_alignment_use_whitebox` (plus optional `pd1_whitebox_tau_days`) so the alignment driver runs the kon/koff/internalisation ODE instead of the legacy filter.
3. Use the updated `scripts.fit_observables` (now supports `--scenarios ...`) to tune the PD‑1 parameters jointly across A-series before moving to true white-box integration.

### Step 3 – Identify the T-cell density driver *(scaffolded; pausing)*

See `docs/tcell_driver_plan.md` for the detailed notes. Summary:
1. Quantify surrogate vs MATLAB trends using the flat-debug CSVs (`t_cells`, `tcell_density_per_ul`).
2. Apply/tune the `tcell_alignment_*` follower knobs already exposed in `alignment_driver_block` to regain basic trend fidelity.
3. When ready, port the minimal `nT1/aT1/T1` ODE slice into the module behind a future `tcell_alignment_use_whitebox` flag.

### Step 4 – Multi-scenario joint fitting *(deferred)*

Once the structure is stable (PD‑1 toggle + clarified T-cell driver):
1. Select fit parameters:
   - Global knobs: geometry weights, PD‑1 kon/koff scales, PK scaling.
   - Optional per-scenario adjustments (e.g. IC tweaks) only if strictly necessary.
2. Define a multi-scenario cost (weighted sum over A1–A6 observables), emphasising PD‑1 accuracy for high-dose regimens.
3. Run `scripts/fit_observables.py` with the multi-scenario objective.
4. Write the resulting parameters back into `parameters/example1_parameters.json` and the snapshot CSV, then rerun:

```
python -m scripts.validate_surrogate \
  --scenarios A1 A2 A3 A4 A5 A6 \
  --ic-mode snapshot \
  --numeric-gates
```

Success target: rel_RMSE for `tumour_volume_l`, `pd1_occupancy`, and `tcell_density_per_ul` < 0.1 (ideally < 0.05) across A-series. At that point A-series is “green” and white-box exporter work becomes incremental rather than firefighting.

---

# White‑box alignment plan

## White‑box alignment plan (PD‑1, tumour, and T cells)**

> Going forward, we will replace the current grey‑box alignment driver with a fully white‑box implementation of the PD‑1 checkpoint and tumour/T‑cell dynamics, while *preserving* the existing frozen‑snapshot infrastructure (equations.txt, JSON parameters, SimBiology‑derived snapshots, and the validate_surrogate tooling). Concretely, we will implement explicit ODEs in Python for (i) antibody PK in the central/peripheral/tumour/LN compartments, (ii) synaptic PD‑1/PD‑L1/PD‑L2 binding and antibody engagement, and (iii) cancer and T‑cell dynamics, including the tumour volume rule and the PD‑1‑dependent killing term, using the published QSP‑IO equations as the primary source of truth. These white‑box modules will live in the alignment layer (e.g. as an extended `alignment_driver_block`) and will read parameters directly from the snapshot/JSON parameter sources. They will *not* change the existing equations.txt or snapshot structure; instead, they will run alongside the frozen model and overwrite only the observables we care about (PD‑1 occupancy, tumour volume, T‑cell density, and derived outputs). This keeps current A‑series/B‑series workflows working, while giving us a transparent, easily extensible Python implementation for future biology changes and Bayesian optimisation.

### Step 1 – Carve out a clean “white‑box alignment” mode

**Goal:** Be able to switch between *snapshot‑only*, *grey‑box*, and *white‑box* drivers without touching the rest of the code.

**Tasks:**

1. **Add a mode switch** in the alignment driver, controlled purely by parameters:

   * e.g. `alignment_mode ∈ {0: snapshot_only, 1: greybox, 2: whitebox}` or two booleans:

     * `pd1_alignment_use_whitebox`
     * `tcell_alignment_use_whitebox`
2. **Keep the public API stable**:

   * `simulate_frozen_model(..., module_blocks=("alignment_driver_block",))` remains the only entry point.
   * The block decides internally (from context/parameters) which branch to run.
3. **Freeze today’s behaviour as mode 0/1**:

   * Mode 0: do nothing, purely snapshot outputs.
   * Mode 1: current grey‑box alignment_driver_block (pd1 filter + geometry/t‑cell followers) – keep it around so you have a working fallback.
4. **Add explicit “white‑box state namespace”**:

   * Reserve keys like `C_hat`, `T_hat`, `H_PD1_hat`, `V_T_hat`, `pd1_YY1_hat`, etc., so you never conflict with frozen snapshot species.
   * All white‑box dynamics write into `*_hat` and then finally map to observable names (`tumour_volume_l`, `pd1_occupancy`, `tcell_density_per_ul`).

This step is mostly plumbing; once it’s in, you can iterate on biology without destabilising A1–A6 runners.

---

### Step 2 – Implement white‑box PD‑1 checkpoint in the alignment driver

**Goal:** Replace the PD‑1 grey‑box filter with an explicit receptor‑binding ODE and Hill‑type inhibition, using the QSP‑IO equations.

#### 2.1 PD‑1 / PD‑L1 / PD‑L2 kinetics

From the QSP‑IO immune checkpoint module:

* Let

  * (Y) = PD‑1 on T cells
  * (Y_1) = PD‑L1 on cancer cells
  * (Y_2) = PD‑L2
  * (A) = anti‑PD‑1 antibody
  * (A_1) = anti‑PD‑L1 antibody

Core synapse reactions:

[
\frac{d[YY_1]}{dt} = k_{\mathrm{on}}^{YY_1}[Y][Y_1] - k_{\mathrm{off}}^{YY_1}[YY_1]
]

[
\frac{d[YY_2]}{dt} = k_{\mathrm{on}}^{YY_2}[Y][Y_2] - k_{\mathrm{off}}^{YY_2}[YY_2]
]

[
\frac{d[YA]}{dt} = \frac{2}{\gamma_T}k_{\mathrm{on}}^{YA}[Y][A] - k_{\mathrm{off}}^{YA}[YA]
]

[
\frac{d[YAY]}{dt} = \frac{\chi}{A_{\mathrm{syn}} d_{\mathrm{syn}} N_A}k_{\mathrm{on}}^{YAY}[YA][Y] - 2k_{\mathrm{off}}^{YAY}[YAY]
]

(and analogous equations for PD‑L1 antibody, if you need them).

**Implementation sketch:**

* Inside `alignment_driver_block` (white‑box mode):

  * Maintain *internal* states: `[YY1_hat]`, `[YY2_hat]`, `[YA_hat]`, `[YAY_hat]`.
  * Take the antibody concentration in tumour, `[A]_T`, from:

    * Either the frozen PK states (`V_T.nivolumab`) if they’re sensible; or
    * Your existing PK surrogate (which already reads the dose schedule and evolves `[A]_C`, `[A]_T`).
  * Use parameter values from your snapshot / JSON:

    * `kon_PD1_PDL1`, `koff_PD1_PDL1`, `kon_PD1_PD1Ab`, `koff_PD1_PD1Ab`, `chi`, `A_syn`, `pd1_synapse_depth_um`, `gamma_T`, `N_A`, etc. (name mapping to QSP‑IO’s (k_{\text{on}}^{YY1}, A_{\text{syn}}, d_{\text{syn}}, \gamma_T) is straightforward).
  * Integrate these ODEs *locally* in the driver using a simple explicit step with the current time step (`dt = t - t_prev`); you don’t need a nested solver, because the step‑to‑step changes over 1 day are small and you are already being updated at every solver call.

#### 2.2 PD‑1 functional inhibition (H_{\mathrm{PD1}})

From the cancer module:

[
H_{\mathrm{PD1}} = \frac{X^n}{X^n + PD1_{50}^n}
]

with

[
X = A_{\mathrm{syn}}[YY_1] \quad \text{(number of PD‑1/PD‑L1 complexes in the synapse).}
]

**Implementation sketch:**

* Define `pd1_whitebox_raw_occ = X = A_syn * YY1_hat`.
* Then compute `H_PD1_hat` via the Hill function above, using `PD1_50` and Hill coefficient `n` from your snapshot.
* Finally set:

  * `context["pd1_occupancy"] = H_PD1_hat` (overriding frozen snapshot value).
  * Log `pd1_whitebox_raw_occ` and `H_PD1_hat` to `--dump-flat-debug` for sanity checks.

Once this is in place, you can see in the flat debug dumps whether the PD‑1 occupancy waveform matches MATLAB qualitatively before you even start global fitting.

---

### Step 3 – Implement white‑box tumour + T‑cell driver for A‑series observables

**Goal:** Generate realistic dynamics for `tumour_volume_l` and `tcell_density_per_ul` in Python, using the published cancer and T‑cell equations, but *without touching* the frozen ODE state vector.

#### 3.1 Cancer cell dynamics (C)

QSP‑IO’s cancer module:

[
\frac{dC}{dt}
= k_{\text{growth}} C \left(1-\frac{C}{C_{\max}}\right)
-\left(k_{\text{innate}} + k_{T_{\text{cell}}} \frac{T}{T+C}(1-H_{\mathrm{PD1}})\right)C
]

* (k_{\text{growth}}): maximal cancer growth rate
* (C_{\max}): carrying capacity
* (k_{\text{innate}}): innate‑immunity death rate
* (k_{T_{\text{cell}}}): T‑cell‑mediated killing rate
* (T): effector T cells in tumour
* (H_{\mathrm{PD1}}): PD‑1 inhibition fraction defined above.

**Implementation sketch:**

* Introduce a white‑box state `C_hat` in the alignment driver.
* Initialize `C_hat` from a sensible value:

  * Either from snapshot `C1` at t=0; or
  * From the target‑volume IC procedure (invert initial tumour volume with your cell volume).
* At every call of `alignment_driver_block` (white‑box mode):

  * Read `H_PD1_hat` from Step 2 and `T_hat` (see below).
  * Integrate `C_hat` forward with a simple ODE step using the equation above.
* Keep all ODE parameters (`k_growth`, `C_max`, `k_innate`, `k_T_cell`) in your existing parameter JSON/snapshot files (many of them already exist there with the same names).

#### 3.2 T‑cell dynamics (minimal version)

You don’t necessarily need the full, multi‑compartment T‑cell system at first. A minimal follower that reacts to PD‑1 and tumour load is often enough for A‑series alignment:

[
\frac{dT}{dt}
= k_{\text{exp}} , T , f_{\text{activation}}(C)

* k_{\text{death}}^T T

* k_{\text{exhaust}} H_{\mathrm{PD1}} T
  ]

* (f_{\text{activation}}(C)) can be a simple saturating function (e.g. (C/(C + C_{50}))) driven by tumour burden.

* All three rates (k_{\text{exp}}, k_{\text{death}}^T, k_{\text{exhaust}}) can live in the parameter JSON as `tcell_alignment_*` knobs (you already have a few T‑cell follower parameters wired through).

If you later want full fidelity, you can port the more detailed activated‑T equations from QSP‑IO, which scale infiltration with tumour volume:

[
\frac{dT}{dt} = q_T^{\text{in}} V_T T_C - k_{\text{death}}^T T + \ldots
]

…but I’d start with the minimal version until A1–A6 look right.

#### 3.3 Tumour volume and T‑cell density observables

From QSP‑IO, tumour volume is:

[
V_T
= V_{\text{cancer}} C_{\text{total}} + V_{T_{\text{cell}}}(T_{\text{total}} + T_{\text{reg}}) + V_{T_{\min}}
]

**Implementation sketch:**

* In the alignment driver (white‑box mode), compute:

  ```text
  V_T_hat = V_cancer * C_hat + V_Tcell * (T_hat + Treg_hat) + V_Tmin
  ```

* Map this into the context:

  ```python
  context["tumour_volume_l"] = V_T_hat
  context["V_T"] = V_T_hat  # if you want downstream algebra to see it
  ```

* For intratumoural T‑cell density (per µL):

  ```text
  tcell_density_hat = T_hat / (V_T_hat * 1e6)
  context["tcell_density_per_ul"] = tcell_density_hat
  ```

* `V_cancer`, `V_Tcell`, `V_Tmin` already have analogues in your parameter files (cell‑volume parameters and minimal tumour volume). If any is missing, add it to JSON and ensure it propagates into snapshot `parameters.csv`.

This gives you full white‑box control over the three observables that currently fail the numeric gates.

---

### Step 4 – Calibrate and lock A‑series with the white‑box layer

**Goal:** Bring A1–A6 tumour volume, PD‑1 occupancy, and T‑cell density within your numeric‑gate thresholds *using only white‑box parameters*, then freeze this as the new baseline.

**Tasks:**

1. **Wire the fitter to white‑box knobs only**:

   * Restrict `scripts/fit_observables.py` to parameters like:

     * `pd1_whitebox_*` (kon/koff effective scalings, PD1_50, Hill n)
     * `pd1_pk_*` if you keep a PK surrogate
     * `tcell_alignment_*`
     * `V_cancer`, `V_Tcell`, `V_Tmin` if needed.
   * Use `--scenarios A1 A3 A5 A6` to force global calibration across dose levels.
2. **Optimisation loop**:

   * For each candidate parameter set, run
     `python -m scripts.validate_surrogate --scenarios A1 A3 A5 A6 --ic-mode snapshot --module-block alignment_driver_block`
     and read the rel_RMSE/maxRE metrics from `metrics.csv`.
   * Minimise a weighted sum of errors on:

     * `tumour_volume_l`, `pd1_occupancy`, `tcell_density_per_ul`.
3. **Check internal signals**:

   * For the final candidate, dump
     `--dump-flat-debug 5` for A1 and A6 and verify:

     * `pd1_whitebox_raw_occ` and `H_PD1_hat` follow the expected MATLAB shape.
     * `C_hat`, `T_hat`, `V_T_hat` roughly resemble the reference tumour/T‑cell trajectories.
4. **Bake parameters into snapshots**:

   * Once you’re happy, copy the tuned values into:

     * `parameters/example1_parameters.json`
     * `artifacts/matlab_frozen_model/example1/parameters.csv`
   * Re‑run `python -m scripts.validate_snapshot ...` and `python -m scripts.run_alignment_suite --scenarios A1-6` to regenerate references and confirm that gates pass without any CLI overrides.

At that point, the A‑series is driven entirely by the **Python white‑box alignment layer**, while the rest of the model (antigen, APC, Tregs, etc.) still comes from the frozen snapshot. That’s exactly the halfway house you wanted: you can now **edit the key biology in Python** (add microenvironment terms, metabolic sinks, Bayesian priors…) without touching the existing equations.txt/SimBiology export machinery, and later you can decide whether to migrate more modules over.

If you’d like, I can next help you turn this into a concrete `TODO` checklist for your team (who owns which step, expected time, and risk points).


### What ODEs and parameters do you actually need?

You already *have* the math in three places:

1. **QSP‑IO tutorial paper (Sové et al. 2020)** – gives closed‑form equations for the cancer module, checkpoint module, antigen module, and PK.
2. **The Jafarnejad NSCLC paper** (the example1 model is “Jafarnejad + minor tweaks”). 
3. **Your MATLAB model** – the same equations appear in `example1.m` + the module functions, and `getequations(model)` has already been exported to `equations.txt` in your repo.

Below are the key bits you’ll want to treat as canonical when you implement the white‑box slice.

---

### 1. Cancer module (tumour cell ODE and PD‑1 inhibition)

From the QSP‑IO tutorial, the tumour cell count (C) in the tumour compartment obeys a logistic growth + immune‑mediated death ODE: 

[
\frac{dC}{dt}
= k_{\text{growth}}, C \Bigl(1 - \frac{C}{C_{\max}}\Bigr)

* \Bigl(k_{\text{innate}} + k_{T_{\text{cell}}} \frac{T}{T + C},\bigl(1 - H_{\mathrm{PD1}}\bigr)\Bigr), C
  ]

where

* (k_{\text{growth}}): cancer proliferation rate (day⁻¹)
* (C_{\max}): carrying capacity (cells)
* (k_{\text{innate}}): baseline innate killing rate (day⁻¹)
* (k_{T_{\text{cell}}}): maximal T‑cell‑mediated kill rate (day⁻¹)
* (T): cytotoxic T cells in the tumour (cells)

The PD‑1 “brake” (H_{\mathrm{PD1}}) is a Hill function of an effective PD‑1 engagement signal (X): 

[
H_{\mathrm{PD1}} = \frac{X^n}{X^n + \mathrm{PD1}_{50}^n}
]

* (X): some measure of PD‑1:ligand complexes on T cells (coming from the checkpoint module)
* (\mathrm{PD1}_{50}): half‑maximal occupancy parameter
* (n): Hill exponent (dimensionless)

**Plan for white‑box:**

* Read (k_{\text{growth}}, C_{\max}, k_{\text{innate}}, k_{T_{\text{cell}}}, \mathrm{PD1}_{50}, n) directly from `parameters/example1_parameters.json` / snapshot `parameters.csv`.
* Implement this ODE literally in Python for the tumour clone used in A‑series (C1).
* Define (X) in terms of synapse‑level species (PD‑1 bound by PD‑L1 / PD‑L2 / drug), which comes from your checkpoint ODEs below.

---

### 2. Immune checkpoint (PD‑1 / PD‑L1 / PD‑L2 + antibodies)

The checkpoint module lives in a synapse compartment and tracks PD‑1, ligands, and antibodies. The tutorial gives explicit ODEs for complexes such as PD‑1:PD‑L1, PD‑1:PD‑L2, PD‑1:drug, and bivalent complexes like drug:PD‑1:PD‑1.

In compact form (using the notation of the paper):

* (Y) = PD‑1
* (Y_1) = PD‑L1, (Y_2) = PD‑L2
* (A) = anti‑PD‑1 antibody (nivolumab)
* (A_1) = anti‑PD‑L1 antibody

Core examples (single‑arm and bivalent drug binding):

[
\frac{d[YA]}{dt}
= \frac{2}{\gamma_T} k_{\text{on}}^{YA},[Y][A]

* k_{\text{off}}^{YA},[YA]
  ]

[
\frac{d[YAY]}{dt}
= \frac{\chi}{A_{\text{syn}} d_{\text{syn}} N_A}
,k_{\text{on}}^{YAY},[YA][Y]

* 2 k_{\text{off}}^{YAY},[YAY]
  ]

(and similarly for PD‑1:PD‑L1, PD‑1:PD‑L2, PD‑L1:drug, PD‑L1:drug:PD‑L1, etc.)

Where:

* (k_{\text{on}}^{\cdot}, k_{\text{off}}^{\cdot}): 2D/3D binding rate constants
* (\gamma_T): tumour interstitial volume fraction
* (\chi): cross‑arm efficiency for bivalent antibodies
* (A_{\text{syn}}): synapse surface area
* (d_{\text{syn}}): synapse gap thickness
* (N_A): Avogadro’s number

**Plan for white‑box:**

* Decide a minimal set of synapse states to carry explicitly (e.g. free PD‑1 (Y), free PD‑L1 (Y_1), PD‑1:PD‑L1 complex ([YY_1]), PD‑1:drug ([YA]), bivalent ([YAY])).
* Translate the corresponding ODEs from the paper/`equations.txt` into a Python function `rhs_checkpoint(t, state, params)` using your unit‑normalised kon/koff from `units.py`.
* Define (X) for the PD‑1 Hill function as some linear combination of PD‑1–inhibiting complexes (e.g. total PD‑1 bound by ligand vs bound by drug), mirroring what Jafarnejad/QSP‑IO does in its repeated assignments.

You don’t have to re‑invent the biology here, just port the equations in a transparent way and let the JSON parameters drive the numbers.

---

### 3. Antibody PK (nivolumab) – 4‑compartment model

The PK module is a standard 4‑compartment antibody model with central (C), peripheral (P), tumour (T), and lymph node (LN) compartments. The tutorial gives the ODEs in terms of compartment volumes (V_\cdot), interstitial volume fractions (\gamma_\cdot), and flows (Q_\cdot):

Representative equations:

[
V_C \frac{d[A]*C}{dt}
= \sum*{i=P,T,\text{LN}} Q_i \Bigl(\frac{[A]_i}{\gamma_i}

* \frac{[A]_C}{\gamma_C}\Bigr)

- Q_{LD}\frac{[A]*{LN}}{\gamma*{LN}}

* k_{cl}[A]_C V_C
  ]

[
V_T \frac{d[A]_T}{dt}
= Q_T \Bigl(\frac{[A]_C}{\gamma_C} - \frac{[A]_T}{\gamma_T}\Bigr)

* Q_{LD}\frac{[A]_T}{\gamma_T}
  ]

(and analogous equations for (A_P), (A_{LN})).

**Plan for white‑box:**

* Implement these four ODEs as a small PK subsystem `rhs_pk(t, pk_state, params)`.
* Drive it with the same dose schedule you already pass into the snapshot engine (you already have `ScheduledDose` etc).
* Feed the tumour PK state ([A]_T) into the checkpoint module as the source for synapse drug concentration.

Because the PK structure is simple and parameters are well defined in your JSON (Vc, Vp, Vt, Qp, Qt, Qln, k_cl, gamma_*), this is one of the easiest modules to white‑box.

---

### 4. Tumour volume & T‑cell density observables

The paper defines tumour volume as a function of cancer cells and cell volume; RECIST diameter is derived assuming a sphere.

You can keep your current geometry logic but make it transparently white‑box:

* Choose a consistent per‑cell volume (v_{\text{cell}}) (or reuse QSP‑IO’s `vol_cell` parameter).

* Define tumour volume (V_T = (C_{\text{live}} + C_{\text{dead}}),v_{\text{cell}}) (plus any dead‑cell swelling factor you already use).

* Compute RECIST diameter

  [
  D_T = 2 \left(\frac{3 V_T}{4 \pi}\right)^{1/3}
  ]

* Define intratumoural T‑cell density as
  [
  \rho_T = \frac{T}{V_T \times 10^6}
  ]
  in “cells per µL”, which matches your current `tcell_density_per_ul` convention.

These are already how you *intend* to compute the observables; white‑boxing just means you compute them directly from the ODE states instead of trying to reconstruct them post‑hoc from a black‑box snapshot.

---

### 5. Where these parameters live in MATLAB / your repo

* In **MATLAB/QSP‑IO**, the parameters are entries in the JSON parameter files, loaded via `load_parameters`, and attached to the SimBiology model as `sbio.Parameter` objects with `Value` and `Units`.
* In your **Python port**, those same values already flow into:

  * `parameters/example1_parameters.json` (canonical author‑side source)
  * the frozen snapshot CSV `artifacts/matlab_frozen_model/example1/parameters.csv`
  * and then into `FrozenModel.parameters` via `load_frozen_model`.

For white‑boxing, you don’t need any *new* parameter source; you just need to:

1. Decide which existing names map to the symbols in the equations above (e.g. `k_C1_growth → k_growth`, `C_max → C_max`, `k_C_T1 → k_Tcell`, `PD1_50 → PD1_50`, `kon_PD1_PDL1` etc.).
2. Pull them out of `model.parameters[...]` in your new white‑box module.
3. Stop introducing separate “alignment_*” knobs except where you explicitly want an extra phenomenological degree of freedom.



## Current status (2025‑11‑13) & next direction

- Steps 0–1 are complete and repeatable (A-series baselines + driver diagnostics in place).
- Steps 2–3 now have the necessary hooks (`pd1_alignment_use_whitebox`, `tcell_alignment_*`) but haven’t been fully calibrated—metrics still show large PD‑1/T-cell errors.
- Rather than iterating further on the grey-box layer, we plan to **pause all incremental tuning and migrate directly to a full white-box implementation** (export complete PD‑1/T-cell subsystems from MATLAB and retire the alignment driver logic).
- Action: keep the above plan as historical context, but start drafting the new “full white-box” roadmap in a fresh section/file.
