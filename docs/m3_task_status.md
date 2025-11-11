# M3 Frozen Snapshot Recovery Checklist

This note tracks the six follow-up items that were proposed for the end of M3.  Each
entry captures the current implementation status and the supporting evidence.

## 1. Re-evaluate aliases after every event and dose

**Status:** ✅ Completed.

`simulate_frozen_model` rebuilds the context through `reconcile`, which now calls
`inject_output_aliases` after replaying repeated assignments, module blocks, and
algebraic rules.  The same helper is used during the sample loop, ensuring that
post-event and post-dose snapshots expose refreshed derived outputs.【F:src/offline/simulation.py†L506-L528】【F:src/offline/simulation.py†L770-L797】

## 2. Sample output columns from the refreshed context

**Status:** ✅ Completed.

When exporting sample rows, the simulator reconstructs each context with module
blocks and alias injection before reading tumour volume, PD-1 occupancy, and the
T-cell density columns.  The exported arrays are therefore sourced directly from
the up-to-date context rather than back-computing from stale state vectors.【F:src/offline/simulation.py†L770-L797】

## 3. Runtime PD-1 bridge module

**Status:** ✅ Completed.

`pd1_bridge_block` is implemented and registered so the bridge can be enabled by
name.  Its logic combines all relevant compartment concentrations with the
`gamma_*` weights and falls back to tumour or central concentrations if no weights
are provided.【F:src/offline/modules/switches.py†L49-L117】

## 4. Tumour geometry refresh module

**Status:** ✅ Completed.

`tumour_geometry_block` back-fills tumour volume when the snapshot carries only cell
counts and then derives tumour diameter and intratumoural T-cell density.  The block
is also registered for runtime use alongside the PD-1 bridge.【F:src/offline/modules/switches.py†L80-L117】

## 5. Quick regression on snapshot-derived observables

**Status:** ❌ Outstanding.

A manual smoke test shows that PD-1 occupancy remains effectively flat when the
bridge is enabled for the `example1` snapshot, and the latest `A1` numeric-gate run
still fails with large discrepancies in tumour volume, PD-1 occupancy, and the
T-cell density columns.  The surrogate/regression deltas therefore have not been
brought back under control.【1647b9†L1-L2】【db9318†L1-L33】

## 6. Regression tests for bridge and geometry dynamics

**Status:** ❌ Outstanding.

The current unit tests exercise only the algebraic combination logic for the PD-1
bridge and the static geometry calculations; there are no automated scenarios that
prove PD-1 occupancy varies after a bolus or that T-cell density tracks changing
volume during integration.  Additional regression cases are still required to lock
in the intended dynamics.【F:tests/test_module_blocks.py†L11-L78】

