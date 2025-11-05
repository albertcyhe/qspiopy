# QSP‑IO Python Runtime — Software Architecture

This document captures the refactored structure of the MATLAB‑free (offline) runtime, its responsibilities, and guidelines for extending the code base.

## Goals

- **Modularity** – isolate concerns (snapshot loading, numerical runtime, domain entities) to make changes local and traceable.
- **Testability** – make it straightforward to unit‑test snapshot parsing, dosing schedulers, and integration logic in isolation.
- **Extensibility** – provide public interfaces that can support new solvers, observables, or scenario definitions without rewriting the core engine.

## High‑level layout

```
src/offline/
├── __init__.py          # Public API (re-exports)
├── entities.py          # Core dataclasses shared across modules
├── snapshot.py          # Snapshot loading and FrozenModel construction
├── simulation.py        # Numerical runtime (dosing, events, integration)
└── frozen_model.py      # Backwards-compatible façade exporting the old API
```

### `entities.py`

Holds the lightweight data structures shared across the runtime:

- `ScenarioResult` – immutable container for simulation outputs (with optional extra observables).
- Domain dataclasses (`SpeciesEntry`, `ReactionEntry`, `EventEntry`, `DoseEntry`, …) used both by the loader and the runtime.
- `ScheduledDose` – sorted dosing entries consumed by the runtime’s scheduler.

This module is intentionally dependency-light (NumPy/Pandas/SymPy only) to keep the types usable from tests without pulling in heavy solver code.

### `snapshot.py`

Responsible for translating a frozen snapshot directory into a `FrozenModel` instance:

- Loads CSV metadata (species, parameters, rules, reactions, events, stoichiometry, …).
- Compiles algebraic/ODE expressions via SymPy (`CompiledExpression` instances).
- Normalises units and derived quantities.
- Provides cached `load_frozen_model(...)`, `snapshot_digest(...)`, and `sha256_file(...)` helpers for provenance tooling.
- The `FrozenModel` class exposes APIs for:
  - Creating the initial state / context.
  - Evaluating rules and reactions.
  - Applying doses and synchronising state vectors.

All filesystem interaction is concentrated here; higher layers treat `FrozenModel` as an immutable runtime description.

### `simulation.py`

Hosts the numerical runtime that turns a `FrozenModel` into trajectories:

- Constants (`SEMANTICS_VERSION`, solver tolerances, event log schema).
- Dosing utilities (`_enumerate_dose_times`, `_build_scheduled_doses`), including injectable custom dosing.
- Diagnostics (`_perform_t0_quick_check`) and sampling helpers.
- `simulate_frozen_model(...)` – the public entry point used by CLI tools/tests. It handles:
  - Solver configuration (tolerances, step sizes, seed banners).
  - Event integration (immediate/delayed queue handling).
  - Output sampling on configurable grids and optional extra observables.

The runtime consumes only the `FrozenModel` API, removing any CSV/JSON reload logic from this layer.

### `frozen_model.py`

Kept as a façade for backwards compatibility. External callers (e.g. CLI scripts) continue to import `simulate_frozen_model` or `EVENT_LOG_FIELDS` from `src.offline.frozen_model` without change. Internally it simply re-exports definitions from the new modules.

## Public API

All upstream tooling should use the re-exported symbols in `src.offline`:

```python
from src.offline import simulate_frozen_model, load_frozen_model, ScenarioResult, DoseEntry
```

This shields callers from future internal rearrangements.

## Testing guidance

- Loader unit tests can target functions in `snapshot.py`, asserting that CSV fixtures produce expected `FrozenModel` metadata and hashes.
- Simulation behaviour (events, dosing, sampling) can be exercised via `simulate_frozen_model` with controlled `custom_doses`, seeds, and `extra_outputs`.
- Because `ScenarioResult` lives in `entities.py`, its `.to_frame()` helper can be validated independently.

## Future extensions

- Add a `solvers/` subpackage if alternative integrators or tolerance strategies are required.
- Introduce dedicated `dosing.py` / `events.py` once scenario variability grows beyond PD‑1/PD‑L1 use cases.
- Provide typed configuration objects (Pydantic models) for CLI scenario definitions to complement the current JSON payloads.
