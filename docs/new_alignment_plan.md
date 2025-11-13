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
