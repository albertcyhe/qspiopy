# Alignment Driver Signal Flow

This note captures the pieces a newcomer needs for handoff step 4.2: how the A-series scenarios move data from frozen SimBiology snapshots through `alignment_driver_block` and into the observables we compare to MATLAB.

## Runtime plumbing (text diagram)

```
MATLAB snapshot (artifacts/matlab_frozen_model/<snapshot>)
    └─ parameters.csv / reactions.csv / state.csv
          ↓ load_frozen_model()
Scenario spec (scripts/scenario_registry.A1)
    └─ doses + context_outputs + module_blocks
          ↓ simulate_frozen_model(..., custom_doses=A1_DOSES)
Segment integrator (src/offline/segment_integrator.py)
    └─ produces state trajectory + module context callbacks
          ↓ alignment_driver_block (src/offline/modules/alignment_driver.py)
                ├─ reads pk inputs: aPD1, V_T.* context
                ├─ updates pd1_filter_* / geom_* fields on context
                └─ writes observables: pd1_occupancy, tumour_volume_l, tcell_density_per_ul
          ↓ ScenarioResult (src/offline/entities.py)
                ├─ samples time grid -> *_surrogate.csv
                └─ exposes raw_contexts for flat-debug capture
```

Key takeaway: all PD‑1 and geometry overrides live inside the module, so the snapshot equations remain untouched while we iterate on surrogate behaviour.

## Instrumentation checklist

1. Invoke validation with context capture (example command):

   ```bash
   python -m scripts.validate_surrogate \
     --scenarios A1 \
     --ic-mode snapshot \
     --module-block alignment_driver_block \
     --numeric-gates \
     --dump-flat-debug 50
   ```

2. The enhanced `_dump_flat_debug` now prints the usual span summary and writes `artifacts/validation/A1_flat_debug_<timestamp>.csv`. Columns include both raw snapshot states (`V_T`, `aPD1`) and alignment driver internals (`pd1_filter_*`, `geom_*`).

3. Plotting overlays live in `artifacts/validation/`; run `scripts/dev_pd1_driver_compare.py` to generate per-scenario PK/binding/occupancy figures straight from the latest flat-debug capture.

### Runtime toggles

- `pd1_alignment_use_whitebox` (plus optional `pd1_whitebox_tau_days`) swaps the legacy filter for the snapshot-derived kon/koff/internalisation ODE.
- `tcell_alignment_*` (`tau_days`, `w_live`, `w_treg`, `offset_cells`, `min_cells`, `occ_supp_coeff`) implements a first-order follower that feeds `tcell_alignment_state` → `tcell_density_per_ul` until the true `nT1/aT1/T1` slice is ported.

These outputs satisfy the “管线已跑通 + 有可视化” requirement from 4.2 and unblock the PD‑1 sandbox work in 4.3.
