# QSP‑IO Python vs MATLAB Alignment — Test Log

## Executed

- **Example 1 / Example 2 regression suite (MATLAB frozen snapshots)**  
  - Runner: `python -m scripts.validate_surrogate --scenarios example1_control example1_treated example2_treated --emit-diagnostics`  
  - artefacts under `artifacts/validation/`: reference/surrogate CSVs, alignment metrics, solver logs.  
  - Result: pass — max relative error across tracked observables ≈7e‑15 (numerical noise), confirming the Python frozen-model runtime remains aligned with the MATLAB-generated snapshot data.

- **A1 (PD‑1 monotherapy, 200 mg Q3W ×4, 12 weeks)**  
  - Runner: `python -m scripts.run_alignment_suite --scenarios A1 --emit-metrics`  
  - artefacts:  
    - Python surrogate: `artifacts/extended_validation/A1_surrogate.csv`  
    - MATLAB replay (reconstituted scripts): `artifacts/extended_validation/A1_reference.csv`  
    - Metrics: `artifacts/extended_validation/alignment_metrics_extended.csv`  
  - Result: fail — relative L2 error on tumour volume ≈1.4e‑1, PD‑1 occupancy max error ≈1.0. Indicates the live MATLAB script replay diverges from the frozen snapshot (likely dose mapping / unit conversion gaps). Needs reconciliation before serving as a comparison baseline.

## In Progress

- **MATLAB replay reconciliation** — align the live MATLAB model execution with the frozen snapshot outputs (example1_control / example1_treated / example2_treated) to within ≤1e‑3 relative L2 on core observables. Current A1 mismatch suggests issues in dosing configuration or state extraction that must be resolved.

## Planned Runs

- **Dose matrix (A2–A6)** — once A1 aligns, extend the runner to the remaining monotherapy regimens and regenerate metrics.
- **Combination therapy (B1–B3)** — introduce PD‑1 + CTLA‑4 dosing after confirming combination species/parameters exist in the frozen snapshot.
- **Population grid (C)** — implement cohort sampler with shared RNG seeds and summarise ORR/DCR deltas vs MATLAB.
- **Sensitivity sweep (D)** — ±20 % one‑at‑a‑time perturbations with Spearman ranking comparison.
- **Steady-state/guardrail checks (E/F)** — zero‑dose stability, high‑dose robustness and pathway knock‑out diagnostics.

These queued scenarios remain pending until the baseline MATLAB ↔ Python parity (≤1e‑3 relative L2 on tracked states) is restored.
