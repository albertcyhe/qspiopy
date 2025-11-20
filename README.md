# QSP Stack — Current State & Next Steps

## Where We Are (M1–M4)
- Hybrid runtime: Python `switches.py` micro‑steps SimBiology `FrozenModel` exports to drive PD‑1, T‑cell, and geometry whiteboxes.
- Works for A1 (single dose), but multi‑dose (A6) remains fragile: SimBiology’s coarse steps (≈0.5 d) collide with PD‑1’s high‑frequency needs, leading to time‑scale mismatch, long pending buckets, and occasional instability even with ramp/chunk tuning.
- Engine is dict/while/side‑effect heavy; no JIT/Autodiff/GPU leverage, so performance and inference are capped.

## Pain Points
- Step‑size mismatch injects “dose shocks” into PD‑1, forcing millions of micro‑steps.
- Chunk/ramp and tolerance tweaks help but do not fully tame A‑series multi‑dose scenarios within reasonable wall‑time.

## What’s Next (M5)
- Rebuild as a **pure Python/JAX QSP engine**:
  - Single time‑integration stack under our control (no SimBiology outer solver).
  - JIT‑friendly, array‑based state and RHS; Autodiff‑ready for calibration/uncertainty.
  - Explicit handling of dosing/events and adaptive stepping without cross‑solver impedance.
- Migrate key models (PD‑1, T‑cell, geometry) into this JAX core, then retire the hybrid bridge.

For background, see `docs/m1_m4_summary.md`. The legacy README remains as-is; use this file as the current orientation point.
