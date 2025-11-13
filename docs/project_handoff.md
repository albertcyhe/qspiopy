# QSP‑IO Snapshot Alignment — Handoff Summary (2025‑11‑13)

This document captures the current state of the SimBiology ↔ Python alignment effort so a new contributor can take over without digging through the whole history.

---

## 1. Context & Goal

- **What we ship**: a Python surrogate that integrates frozen SimBiology snapshots (`artifacts/matlab_frozen_model/<snapshot>`) and reproduces MATLAB reference trajectories (`artifacts/validation/<scenario>_reference.csv`).
- **Why it hurts**: the base solver (M1–M3) is stable and unit‑consistent, but A‑series scenarios still fail the numeric gates because PD‑1 occupancy, tumour volume, and t‑cell density stay “flat” in Python while MATLAB shows slow dynamics.
- **Strategy**: keep the snapshot pipeline untouched (“white box”), add a thin “alignment layer” (runtime modules) only where MATLAB semantics are missing, and gradually retire the grey‑box logic as soon as we can export the true equations.

---

## 2. What’s done (kept and tested)

| Area | Status | Notes |
| --- | --- | --- |
| **Integrator & events (M1)** | ✅ | Segmented integration with ε‑bump, trigger specs, pre/post records (`src/offline/segment_integrator.py`, `simulate_frozen_model`). |
| **Units & parameters (M2)** | ✅ | `src/offline/units.py` normalises everything to day/L/mol/M; dose application goes through `normalise_dose_to_species`. |
| **Snapshot runtime (M3)** | ✅ | `simulate_frozen_model` handles snapshot/target‑volume ICs, event de-bounce, warm starts, CLI toggles; module system in `src/offline/modules/switches.py`. |
| **Instrumentation** | ✅ | `--dump-flat-debug` now captures module internals (PD‑1 filter, geometry stats). New script `scripts/dump_snapshot_semantics.py` mirrors MATLAB’s `show_alignment_drivers.m`. |
| **Grey-box alignment layer** | 🚧 | `alignment_driver_block` introduces a minimal PK/PD‑1/volume ODE, but parameters are placeholder values and waveforms still miss MATLAB shapes. |

Artifacts to know:
- `artifacts/show_alignment_example1.txt`: MATLAB report via `matlab/scripts/show_alignment_drivers.m`.
- `artifacts/show_snapshot_example1_filtered.txt`: matching snapshot dump (filtered for PD‑1/geometry keywords).

---

## 3. Known gaps (must fix before gates turn green)

1. **Exporter parity**: MATLAB exporter currently emits only a subset of reactions/rules in its summary. We need to ensure `reactions.csv`/`equations.txt` contain **every** reaction (including PD‑1 transport/binding) so `FrozenModel.rhs()` stays authoritative.
2. **Runtime overrides**: historical modules (`pd1_bridge_block`, `tumour_geometry_*`) disable the very rules we just exported. For A‑series we should run with **only** `alignment_driver_block` until white-box semantics are restored.
3. **Alignment driver calibration**: the new ODEs (PD‑1 binding + logistic tumour volume) are wired up but unfit—occupancy jumps to 0.13 instantly, volume reacts on day 1, and `tcell_density_per_ul` is off by ~1e6×. We need to fit the explicit parameters (`pd1_occ_*`, `geom_*`) against MATLAB references.
4. **Snapshot parameter coverage**: MATLAB `parameters.csv` now includes the alignment knobs (e.g. `pd1_occ_kon_scale`, `geom_growth_per_day`). Once the white-box exporter is fixed, these should move back into SimBiology so Python no longer hardcodes defaults.

---

## 4. Immediate next steps for a new contributor

1. **Exporter audit** (MATLAB side)  
   - Locate the snapshot exporter (`matlab/scripts/export_matlab_snapshot.m` or equivalent).  
   - Ensure it iterates `model.Reactions` / `model.Rules` / `model.Events` and dumps each entry to CSV/TXT.  
   - Re-export `example1` (and any other scenarios) into `artifacts/matlab_frozen_model/<name>`.  
   - Run `python -m scripts.validate_snapshot artifacts/matlab_frozen_model/example1` and `scripts/dump_snapshot_semantics.py …` to confirm parity with MATLAB outputs.

2. **Alignment driver (white-box form)**  
   - Keep `alignment_driver_block` as the sole PD‑1/geometry module for A1.  
   - Verify `--dump-flat-debug` shows the ODE states (`pd1_alignment_*`) moving smoothly once the exporter is fixed.  
   - Use `scripts/fit_observables.py` to fit PD‑1 parameters first (occupancy L2), then tumour volume / t‑cell density. Write tuned values back to `parameters/*.json` + snapshot CSV.

3. **Re-test numeric gates**  
   - Command: `python -m scripts.validate_surrogate --scenarios A1 --ic-mode snapshot --module-block alignment_driver_block --numeric-gates --dump-flat-debug 20`.  
   - Once A1 passes, exercise A2–A6 and B scenarios; add scenario-specific overrides only if absolutely necessary.

4. **Documentation updates**  
   - Keep this file and `docs/new_alignment_plan.md` in sync with each milestone.  
   - Record new MATLAB exports (keep the `artifacts/show_alignment_*.txt` logs) so diffing stays easy.

---

## 5. Repo map for handoff

| Path | Purpose |
| --- | --- |
| `src/offline/` | Snapshot loader (`snapshot.py`), solver (`simulation.py`), runtime modules (`modules/switches.py`). |
| `scripts/validate_surrogate.py` | CLI entrypoint for scenarios; registers A/B series via `scripts/scenario_registry.py`. |
| `scripts/fit_observables.py` | Grey/white-box parameter fitter (least-squares). |
| `scripts/dump_snapshot_semantics.py` | Snapshot-side dump for PD‑1/geometry semantics. |
| `docs/new_alignment_plan.md` | Living summary of milestones/status (now trimmed; see below). |
| `docs/m4_geometry_status.md` | Historical notes on geometry calibration; keep for reference. |
| `artifacts/matlab_frozen_model/*` | Frozen snapshots (input). |
| `artifacts/validation/*` | MATLAB reference CSVs, Python surrogate outputs, diagnostics. |

---

## 6. Contact checklist

When syncing with MATLAB-side owners, align on:
- Exporter coverage (full reaction/rule dump).  
- Official PD‑1 geometry constants (`geometry_*`, `pd1_occ_*`).  
- Whether PD‑1 occupancy requires additional MATLAB-only modules (if so, expose their math so we can reproduce it in Python).

Once these pieces are in place, the remaining work is mainly calibration + regression tests. Good luck!
