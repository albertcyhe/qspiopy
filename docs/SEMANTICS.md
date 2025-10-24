# Frozen SimBiology Execution Semantics

This document captures the Python runtime rules that mirror MATLAB SimBiology for the frozen tutorial snapshots.

## Snapshot Contract

* `equations.txt` is the authoritative source for ODEs and Fluxes. No post-processing beyond token sanitisation is permitted.
* Metadata exported alongside each snapshot – `configset.json`, `species.csv`, `rules.csv`, `events.csv`, `doses.csv`, `variants.csv` – is validated on load and drives runtime behaviour (units, solver tolerances, scheduling, etc.).
* Reference trajectories record the SHA-256 of `equations.txt` and `configset.json`; replacements must bump the versioned artefacts intentionally.

## Initialisation Order

1. Load compartments, parameters, species, constants, config, and rules.
2. Apply initial assignment rules using the frozen parameter/state context.
3. Apply t = 0 doses prior to starting the integrator.
4. Evaluate repeated assignments and algebraic rules; synchronise the state vector to honour nonnegative/boundary constraints.
5. Run the t = 0 diagnostic – fluxes, ODE RHS, and trigger predicates must match the frozen MATLAB snapshot within the configured tolerance.

## RHS Evaluation Pipeline

Per RHS evaluation (`FrozenModel.rhs`):

1. Build the runtime context from solver state (`y`) and static parameters/compartments.
2. Evaluate the repeated-assignment graph in topological order until convergence.
3. Solve algebraic rules (symbolic lambdify where available, otherwise assert residuals within tolerance) without mutating constants.
4. Cache reaction flux values and expose them through the context.
5. Assemble the derivative vector from `ODEs:` entries and overlay any rate rules.
6. Enforce `BoundaryCondition` and `NonNegative` flags in the resulting derivative/state projections.

## Event & Dose Scheduling

* Events keep the original SimBiology order (`event_index`) for same-time triggers; delays queue new `ScheduledEvent` instances for future processing.
* Delayed events, as well as bolus doses at triggering instants, are applied before the solver restarts at `t_event + ε`.
* Doses emitted from snapshots take precedence; a fallback nivolumab schedule is only used when the snapshot is silent and the scenario requests `anti_pd1` therapy.
* All state changes immediately re-run repeated assignments and algebraic rules before re-entering the integrator.

### Event Log Schema

Diagnostics mode writes an event log per scenario (`events_<scenario>_python.csv`). Each row conforms to the schema below and captures both immediate and delayed executions:

| Column | Meaning |
| --- | --- |
| `scenario` | Scenario identifier supplied by the executor |
| `event_index` | Original SimBiology event index (`index_in_model`) |
| `time_fire` | Simulation time at which assignments executed |
| `time_trigger` | Time the trigger condition first evaluated true |
| `delay` | Evaluated delay (same units as `time_fire`) |
| `type` | `"immediate"` or `"delayed"` |
| `assignments` | Canonicalised assignment list `Target=Expr` separated by `;` |

## Constraints & Diagnostics

* `BoundaryCondition` species ignore reaction-driven updates; direct assignments/rate rules remain honoured.
* `NonNegative` species are clamped to zero after each context synchronisation.
* Violations and diagnostic values at t = 0 are exported to `equations_eval_t0.csv` (best effort when the filesystem is read-only).

## CLI & Regression Hooks

* `python -m scripts.validate_surrogate --replicates 1 --max-rel-err 1e-6` is the CI gate; it aborts on t = 0 mismatches or excess fractional error.
* Reference artefacts live in `artifacts/validation` and now include snapshot digests for provenance tracking.
* Unit tests cover trigger-direction parsing, jitter tolerance, and non-linear algebraic rule resolution; extend these to exercise future edge cases (concurrent events, constraint violations, etc.).
