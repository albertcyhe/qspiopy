# QSP‑IO Python vs MATLAB Alignment — Test Log

## Executed

- **Example 1 / Example 2 regression suite (MATLAB frozen snapshots)**  
  - Runner: `python -m scripts.validate_surrogate --scenarios example1_control example1_treated example2_treated --emit-diagnostics`  
  - artefacts under `artifacts/validation/`: reference/surrogate CSVs, alignment metrics, solver logs.  
  - Result: pass — max relative error across tracked observables ≈7e‑15 (numerical noise), confirming the Python frozen-model runtime remains aligned with the MATLAB-generated snapshot data.
- **A1 (registry-driven dry run, Python only)**  
  - Runner: `python -m scripts.run_alignment_suite --scenarios A1 --skip-matlab`  
  - artefacts: `artifacts/extended_validation/A1_surrogate.csv`, header manifest, provenance metadata.  
  - Result: Python candidate generated correctly with new registry wiring; MATLAB reference still pending due to outstanding alignment gaps (see In progress).

- **A1 (PD‑1 monotherapy, 200 mg Q3W ×4, 12 weeks)**  
  - Runner: `python -m scripts.run_alignment_suite --scenarios A1 --emit-metrics`  
  - artefacts:  
    - Python surrogate: `artifacts/extended_validation/A1_surrogate.csv`  
    - MATLAB replay (reconstituted scripts): `artifacts/extended_validation/A1_reference.csv`  
    - Metrics: `artifacts/extended_validation/alignment_metrics_extended.csv`  
  - Result: fail — relative L2 error on tumour volume ≈1.4e‑1, PD‑1 occupancy max error ≈1.0. Indicates the live MATLAB script replay diverges from the frozen snapshot (likely dose mapping / unit conversion gaps). Needs reconciliation before serving as a comparison baseline.

## In Progress

- **MATLAB replay reconciliation (A1)** — align the live MATLAB model execution with the frozen snapshot outputs (example1) to within ≤1e‑3 relative L2 on core observables. Current run still shows large residuals (PD‑1 occupancy and cell counts), pointing to dosing/unit semantics that must be corrected before progressing to the rest of the A-series.
- **Scenario orchestration** — the new registry (`scripts/scenario_registry.py`) and CLI wiring are in place; pending tasks include MATLAB reference regeneration for A1 and subsequent monotonicity/BH checks once numerical agreement is reached.

## Planned Runs

- **Dose matrix (A2–A6)** — once A1 aligns, extend the runner to the remaining monotherapy regimens and regenerate metrics.
- **Combination therapy (B1–B3)** — introduce PD‑1 + CTLA‑4 dosing after confirming combination species/parameters exist in the frozen snapshot.
- **Population grid (C)** — implement cohort sampler with shared RNG seeds and summarise ORR/DCR deltas vs MATLAB.
- **Sensitivity sweep (D)** — ±20 % one‑at‑a‑time perturbations with Spearman ranking comparison.
- **Steady-state/guardrail checks (E/F)** — zero‑dose stability, high‑dose robustness and pathway knock‑out diagnostics.

These queued scenarios remain pending until the baseline MATLAB ↔ Python parity (≤1e‑3 relative L2 on tracked states) is restored.
