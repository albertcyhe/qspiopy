# QSP‑IO Alignment & Tuning Roadmap (A-Series → B-Series)

This document captures the staged plan for bringing the extended scenario matrix
to production quality while keeping the core offline runtime untouched. It
covers the verification loops, acceptance gates, and engineering tasks that
need to be wired into the CLI/CI toolchain.

---

## Verification loops

1. **Semantic alignment (engine ↔ SimBiology)**
   - `t=0` quick check → single low-dose micro scenario → **A1** → **A2–A6** →
     **B-series**. Each stage must pass the numerical thresholds before moving
     forward.
2. **Biological consistency (BH rules)**
   - After every scenario batch, immediately evaluate monotonicity,
     occupancy/peak ranges, and steady-state return. Failures halt progression.

---

## Phase breakdown

### Phase 0 — Tooling prep (one-time)

**Goal:** codify health checks and provenance in the CLI/CI entry points.

- Ensure `simulation.py` exposes the existing P0/P1 self-checks via CLI flags.
- Extend `ScenarioResult` with `semantics_version`, `provenance` (snapshot SHA,
  unit flags, solver config), contract `header`, and `to_frame(order="contract")`
  / `save_csv()` helpers.
- Introduce a lightweight scenario registry (`scripts/scenario_registry.py`)
  describing A/B dosing grids and sampling resolution.

### Phase 1 — A1 baseline (“golden anchor”)

**Goal:** lock down the reference A1 run.

- Numerical acceptance: `rel_L2 ≤ 1e-3`, `maxRE ≤ 5e-3` on tumour/drug/occupancy,
  event time offset ≤ 0.5 h, `t=0` + boundary/non-negative checks clean.
- Use snapshot `TimeUnits` exactly, sample at ≤12 h cadence.
- Centralise mg→amount/concentration conversion in the snapshot layer.
- Perform a tolerance/max-step shrink test (factor 0.5) and ensure drift < 0.5 %.
- Automate the gate via `python -m scripts.validate_surrogate --scenarios A1 --dump-t0 --numeric-gates --emit-metrics`
  (fails fast if any threshold or event timing residual exceeds the limits).

### Phase 2 — A2–A6 matrix

**Goal:** broaden coverage, enforce monotonic exposure-effect behaviour.

- Numerical thresholds identical to A1.
- Biological gates (BH-2/3/4):
  - AUC/Cmax non-decreasing with higher dose/frequency.
  - TGI/BOR ordering matches exposure (allow mild saturation).
  - Occupancy rises during treatment windows and falls after.
- Engineering:
  - Generate custom doses through the shared registry.
  - Extend `validate_surrogate` with `--check-monotonicity` (AUC/TGI/BOR sorting).
  - Emit AUC vs TGI/BOR scatter plots under `artifacts/extended_validation/`.

### Phase 3 — B-series (PD-1 + CTLA-4)

**Goal:** introduce combination therapy without core engine changes.

- Precondition: snapshots include CTLA-4 PK/PD state variables; otherwise
  refresh/export snapshots before tuning.
- Numerical gates: same as Phase 1.
- Biological gates (BH-7/6):
  - Combination TGI ≥ monotherapy, B3 deepest/quickest suppression.
  - IL-6 / IFN-γ peaks increase but remain within physiology-informed bounds.
- Engineering:
  - Registry provides parallel PD-1 + CTLA-4 doses (distinct targets).
  - Extra outputs for dual AUCs, cytokine peaks; CLI `--bh-checks` verifies BH-6/7.
  - Produce side-by-side plots (monotherapy vs combo) and event residual tables.

---

## Metrics & CI gates

- **Per-scenario numeric alignment:** rel_L2 ≤ 1e-3, maxRE ≤ 5e-3 (tumour volume,
  plasma drug, occupancy), event |Δt| ≤ 0.5 h. Failures raise `AlignmentFail`.
- **BIological consistency:**
  - A-series: AUC/TGI/BOR ordering, early BOR advantage for loading dose (A5).
  - B-series: TGI superiority, cytokine peaks within configured limits.
- **CI jobs:**
  1. `validate_surrogate --scenarios A1 --strict` (mandatory for PRs).
  2. `validate_surrogate --scenarios A1,A2,A3,B1 --check-monotonicity --bh-checks`
     (nightly/mainline).

---

## Engineering checklist

1. Add scenario registry (doses, sampling, extra outputs for A/B).
2. Wire registry into `run_alignment_suite` and CLI flags for monotonic/BH checks.
3. Persist contract headers (`scenario*.header.txt`) and use them for alignment.
4. Generate diagnostic plots/reports under `artifacts/extended_validation/`.
5. Emit `dose_audit.csv` from the suite to track target, units, and Δvalues per dose.
6. Optionally add a smoke CI workflow (ruff + mypy + pytest + alignment smoke).
7. Record reference SHA changes in provenance when updating MATLAB exports.

---

## Troubleshooting tips

1. **`t=0` quick check fails:** examine default species dimensions, unit conversion,
   and boundary-condition flags.
2. **Large 24–48 h drift:** re-check dosing target/units and time unit mismatch.
3. **Constant factor offset:** indicates unit scaling error (mass vs concentration).
4. **Event timing errors:** verify trigger direction, delay handling, event index.
5. **Combination mismatch:** confirm CTLA-4 states exist and doses map to the
   correct target species.

---

This roadmap is the source of truth for all A/B tuning work. Update it whenever
new gates or scenarios are added. A change is only “done” when the relevant
CI gates are green and the BH checks pass.

---

## Latest Alignment Gate (A1)

- **Command:** `python -m scripts.validate_surrogate --scenarios A1 --dump-t0 --numeric-gates --output artifacts/validation`
- **Status:** ❌ (numeric gate failed)
- **Violations (`artifacts/validation/A1_surrogate.csv` vs `A1_reference.csv`):**
  - `tumour_volume_l`: rel_L2 = 1.417e-1, maxRE = 1.916e-1
  - `pd1_occupancy`: rel_L2 = 1.0, maxRE = 1.0 (surrogate stays ~0 while MATLAB shows expected dynamics)
  - `tcell_density_per_ul`: rel_L2 = 7.862e-1, maxRE = 7.999e-1
- **Dose audit parity:** numerical deltas now agree, but the surrogate still labels the PD-1 central species as `interpreted_dimension=amount` while MATLAB reports `molarity`, so the audits differ textually (`artifacts/validation/A1_dose_audit.csv` vs `_reference_dose_audit.csv`).
- **Next steps:** inspect tumour module rules (volume ↔ diameter mapping), PD-1 occupancy/repeated assignments, and the T-cell density calculation to close the remaining semantic gaps before expanding to A2–A6/B*.
