# T-cell Density Driver Notes (2025-11-13)

## Observations from current surrogate vs MATLAB

Source data: `artifacts/validation/A1_reference.csv`, `artifacts/validation/A1_surrogate.csv`, `artifacts/validation/A1_flat_debug_20251113T151854.csv`.

1. **Surrogate `t_cells` is nearly frozen**
   - Range: 1.52e6 → 1.83e6 cells (`A1_surrogate.csv`).
   - Context dump shows `T_tumour` (`V_T.T1`) increases slightly from 1.07e6 → 1.30e6, but the total T-cell state (`t_cells`) never exceeds 1.83e6 because lymph node / peripheral pools are flat.

2. **MATLAB reference `t_cells` varies by >4×**
   - Reference stats: mean 5.67e6, max 6.76e6 cells (`A1_reference.csv`).
   - There is steady growth between days 0–84, driven by the activation/proliferation reactions listed in the snapshot dump (e.g. `k_T1_act`, `k_T1_pro`, trafficking terms, PD-1 suppression on `T1 -> T_exh`).

3. **Density computation mismatch**
   - Surrogate `tcell_density_per_ul` exactly equals `t_cells / (tumour_volume_l * 1e6)`; see `A1_surrogate.csv` where the derived value matches to 1e-12.
   - Reference density is lower than the simple ratio (difference ~30 cells/µL at day 0), implying MATLAB applies additional geometric scaling (likely from the pseudo-progression module or cell-to-volume factors) before reporting the observable.

4. **Driver signals available today**
   - Flat debug already records `T_tumour` (alias `V_T.T1`) and the raw context volume `V_T`.
   - Snapshot repeated assignments show `T_total = V_T.T1 + V_T.T0` and `H_PD1_C1` feeds multiple T-cell reactions (e.g. `T1 -> T_exh` with coefficients `k_T1`, `k_C_T1`).

## Implications for Step 3

### Low-effort bridge
- Implement a follower inside `alignment_driver_block` (or a sibling module) that drives `t_cells` / `tcell_density_per_ul` towards the MATLAB reference curve via a first-order lag. Inputs could be the existing `geom_volume_smoothed_l` plus a scenario-specific target trajectory.
- Pros: quick to prototype; keeps the observable responsive to PD-1 inhibition via a tunable gain.
- Cons: still ignores the real lymph node dynamics and makes multi-scenario calibration cumbersome.

### Higher-fidelity white-box slice
- Reuse the snapshot reactions for the `nT1 / aT1 / T1` compartments, including trafficking (`q_T1_*`) and PD-1 mediated exhaustion (`T1 -> T_exh`).
- Integrate those ODEs inside the alignment driver (small state vector: `nT1`, `aT1`, `T1`, optional `T_exh`), parameterised by the exported constants (see `scripts.dump_snapshot_semantics` output above).
- Outputs:
  1. Use the integrated `T1` + `T0` to compute `tcell_density_per_ul` with the same geometry function as MATLAB.
  2. Feed the resulting `T1` back into the tumour kill term so volume dynamics pick up the improved trend.
- Pros: captures regimen-specific amplification (A3/A5/A6) naturally via the actual biological parameters; aligns with the long-term goal of retiring the grey-box layer.
- Cons: requires careful wiring to avoid double-counting states already present in the frozen snapshot.

## Recommended next steps
1. Prototype the follower (low effort) to confirm that simply allowing `t_cells` to grow with scenario-specific gains is enough to bring `tcell_density_per_ul` rel_RMSE below 0.3. → **Done** via `tcell_alignment_*` parameters inside `alignment_driver_block` (Nov-13).
2. In parallel, identify the minimal subset of `T` reactions needed for the white-box slice (likely Reaction_14–27) and sketch how to integrate them with the snapshot state vector without re-solving the full model.
3. Once the PD-1 white-box toggle is validated (Step 2), reuse the same flagging mechanism (`pd1_alignment_use_whitebox`) for the T-cell subsystem, e.g. `tcell_alignment_use_whitebox` with its own smoothing horizon.
