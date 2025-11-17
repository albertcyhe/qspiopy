# PD‑1 Synapse Handoff

This note is the single onboarding document for anyone picking up the PD‑1 synapse alignment work. It summarises what is already done, how to reproduce the current results (Python + MATLAB), and what remains to turn the white‑box module into a production‑ready replacement for the legacy driver.

---

## 1. Where we are today

| Area | Status | Details |
| --- | --- | --- |
| **ODE + solver** | ✅ | `src/offline/modules/pd1_whitebox.py` now implements Reactions 89‑92 exactly as in `artifacts/matlab_frozen_model/example1/equations.txt`. States stay in molecules·µm⁻², fluxes are divided by `syn_T1_C1` area, and integration is delegated to the shared stiff helper (`src/offline/stiff_ode.py`) using the same `SolverConfig` as the frozen model. Unit test: `tests/test_pd1_params.py` evaluates Reaction 89 numerically at t=0 and matches the MATLAB flux. |
| **Parameter loading** | ✅ | `src/offline/modules/pd1_params.py` reads `parameters/example1_parameters.json`, applies the same derivations MATLAB uses (seconds→days, divide by synapse depth, etc.), and returns a structured `PD1Params`. Both the runtime (`alignment_driver_block`) and fitter (`scripts/fit_pd1_whitebox.py`) consume this object. |
| **Diagnostic harness** | ✅ | `scripts/dev_pd1_training_diff.py` compares the Python white‑box with any scenario in `artifacts/training/pd1_whitebox_training.parquet`, writing CSVs under `artifacts/dev/` and printing RMSE for `H_PD1` + block fraction. Example output (current baseline): `pd1_train_0004 H_RMSE=0.403 block_RMSE=0.0423`, `pd1_train_0582 H_RMSE=0.1398 block_RMSE=0.0190`. |
| **MATLAB probing** | ✅ | `matlab/scripts/dev_pd1_probe.m` runs `example1.m`, prints the synapse surface density inferred from the parameters (`PD1 ≈ 396.5358`, `PDL1 ≈ 1783.1862` molecules·µm⁻²), and evaluates Reaction 89. These numbers match what the Python loader now uses. |
| **Overall accuracy** | 🔄 | Block fraction is already within O(10⁻²), but `H_PD1` RMSE is still ≈0.4 (moderate dose) and ≈0.14 (high dose). The remaining gap is due to the Hill transform not matching the MATLAB training curves; this is the focus of the next iteration. |

---

## 2. Assets and file layout

| Path | Purpose |
| --- | --- |
| `artifacts/training/pd1_whitebox_training.parquet` | 600 scenario export from MATLAB (`matlab/scripts/export_pd1_training_suite.m`). Columns: `time_days`, `drug_tumor_molar`, synapse states, and MATLAB’s `pd1_inhibition`. |
| `src/offline/modules/pd1_whitebox.py` | Current PD‑1 module (white‑box) that mirrors SimBiology’s reactions and uses the shared stiff solver helper. |
| `src/offline/stiff_ode.py` | Reusable wrapper around `solve_ivp` (`solve_stiff_ivp` + `integrate_local_system`) shared by the PD‑1 module and the main segmented integrator. |
| `scripts/dev_pd1_training_diff.py` | Lightweight harness to compare the white‑box against any training scenario. |
| `matlab/scripts/dev_pd1_probe.m` | MATLAB probe to print synapse area, PD‑1 densities, and Reaction 89 flux. Run via `/Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab -batch "cd('/Volumes/AlbertSSD/Program/new/qspiopy'); run('matlab/scripts/dev_pd1_probe.m');"`. |
| `scripts/fit_pd1_whitebox.py` | Fitter that ingests the parquet, instantiates `PD1WhiteboxModel`, and optimises selected parameters (`kon`/`koff` scales, PD1_50, internalisation). Uses the same `PD1Params` struct as the runtime. |

---

## 3. Reproducing today’s numbers

1. **Python diff harness**  
   ```bash
   python scripts/dev_pd1_training_diff.py \
       artifacts/training/pd1_whitebox_training.parquet \
       --scenarios pd1_train_0004 pd1_train_0582
   ```
   This writes `artifacts/dev/pd1_compare_<scenario>.csv` with columns:
   ```
   time_days, H_PD1_matlab, H_PD1_python, block_matlab, block_python, H_diff, block_diff
   ```
   and prints the RMSE summary.

2. **MATLAB snapshot probe**  
   ```bash
   /Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab -batch \
       "cd('/Volumes/AlbertSSD/Program/new/qspiopy'); run('matlab/scripts/dev_pd1_probe.m');"
   ```
   Expected output:
   ```
   synapse area (µm^2): 37.8
   PD1 density: 396.535831
   PDL1 density: 1783.186226
   PD1_PDL1 density: 0
   ...
   kon_PD1_PDL1: 0.0583333 (1/(density*day))
   koff_PD1_PDL1: 1.435 (1/day)
   Reaction_89 flux at t0: ~1.56e+06 (molecules/day)
   ```
   Use this to confirm the Python loader continues to align with the snapshot.

3. **Unit tests + lint**  
   ```
   pytest tests/test_pd1_params.py tests/test_stiff_ode_helper.py tests/test_module_blocks.py
   ```
   These cover the Reaction 89 parity test, stiff solver sanity checks, and module plumbing.

---

## 4. What still needs doing

| Priority | Task | Notes |
| --- | --- | --- |
| ⭐️ | **Bring `H_PD1` RMSE below 1e‑2 for pd1_train_0004 + pd1_train_0582** | Use `scripts/dev_pd1_training_diff.py` as the tight loop. Options: adjust `PD1_50`/`n_PD1` via the fitter, verify whether MATLAB’s training curves include extra filtering, or reconcile the `pd1_inhibition` definition if it’s not strictly the Hill output. Document changes once the RMSE drops. |
| ⭐️ | **Re‑run the PD‑1 fitter once the local diff is green** | Command: `python scripts/fit_pd1_whitebox.py --training-path artifacts/training/pd1_whitebox_training.parquet --parameter-file parameters/example1_parameters.json`. Update `parameters/example1_parameters.json` with the tuned values, and keep a copy of the diff CSVs for regression evidence. |
| ⭐️ | **Instrument A‑series runs** | After the training RMSE target is met, re‑enable `python -m scripts.validate_surrogate --scenarios A1 ... --emit-diagnostics ...` and `scripts/dev_pd1_driver_compare.py`. Use the new solver logging (from `stiff_ode.py`) to ensure no segment stalls. |
| ⚪️ | **Documentation clean‑up** | Once PD‑1 white‑box is green end‑to‑end, update `docs/new_alignment_plan.md` to mark the PD‑1 milestone as done, and note the new helper (`stiff_ode.py`) so future modules reuse it. |

---

## 5. Suggested workflow for the next engineer

1. Clone the repo (or sync to the latest branch), ensure you can run `pytest ...` and the diff harness without modification.
2. Iterate on PD‑1 parameters / Hill settings until `scripts/dev_pd1_training_diff.py ...` reports `H_RMSE < 1e-2` for both the moderate and high dose scenarios. Commit the CSVs in `artifacts/dev/` as evidence each time you hit a milestone.
3. Re-run the fitter and update `parameters/example1_parameters.json`.
4. Once the single-scenario RMSE target is met, re-run the A1–A6 diagnostics with `--emit-diagnostics` and capture logs showing the solver is healthy.
5. Only after PD‑1 is stable should we flip the default runtime away from the legacy driver.

If you need to modify MATLAB again, reuse `matlab/scripts/dev_pd1_probe.m` or `matlab/scripts/export_pd1_training_suite.m` so the Python side always has ground-truth data to diff against.

Welcome aboard! The heavy lifting (unit alignment + solver refactor) is already finished; the remaining work is purely about bringing the training and validation RMSE into spec.***

---

## Appendix: legacy Step‑by‑Step plan (for context)

We used the following phased plan earlier in the project. The table below is preserved for historical context; the “Current status” column tells you whether the step is already complete or folded into the new workflow.

| Step | Summary | Current status |
| --- | --- | --- |
| **Step 0 — Orientation** | Review `docs/new_alignment_plan.md`, `equations.txt`, stiff solver code, and the old PD‑1 fitter. | ✅ Covered in Section 2 of this doc. |
| **Step 1 — Extract reusable stiff helper** | Factor out the BDF wrapper (`integrate_local_system`) and add unit tests. | ✅ `src/offline/stiff_ode.py` + `tests/test_stiff_ode_helper.py`. |
| **Step 2 — Refactor PD1WhiteboxModel** | Use the shared helper, implement Reactions 89–92 exactly, remove ad‑hoc clamps, and pass `SolverConfig` down from `alignment_driver_block`. | ✅ See Section 1 and `src/offline/modules/pd1_whitebox.py`. |
| **Step 3 — Fix PD‑1 parameter units** | Align the loader with the snapshot units (`k/kd/Chi`), convert to densities, and feed a structured parameter object into both runtime and fitter. | ✅ Section 1 + `src/offline/modules/pd1_params.py`. |
| **Step 4 — Regression + validation loop** | Run the fitter, inspect training diffs, then re-run A1–A6 with diagnostics. | 🔄 This is the active work tracked in Sections 3–5. |

Keep this appendix as a quick reminder of the original intent; the actionable instructions live in the main sections above.
