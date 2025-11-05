# Changelog

## v0.9.1 (unreleased)

- Added event logging instrumentation (immediate/delayed trigger metadata, assignments) and surfaced the logs through `simulate_frozen_model`, `simulate_reference`, and `scripts.validate_surrogate`.
- Validation now compares the frozen Python execution (`simulate_frozen_model`) directly against the MATLAB snapshot trajectories.
- Introduced scenario-aware diagnostics (`--scenarios`, `--seed`, event log CSVs) and enhanced validation logging with solver banners and worst-error pinpointing.
- Published aggregation utilities (`scripts/summarize_equivalence.py`) generating `alignment_metrics`, `alignment_summary`, `model_scale`, and `event_logs` tables for manuscript figures.
- Documented the event log schema in `docs/SEMANTICS.md` and refreshed README instructions for running diagnostics plus summary scripts.
