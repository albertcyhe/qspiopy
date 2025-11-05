QSP‑IO (Python) — MATLAB SimBiology semantics reproduced in pure Python
=======================================================================

This repository provides a MATLAB‑free execution path for the QSP‑IO immuno‑oncology models by freezing SimBiology projects to disk and replaying their semantics in Python. The goal is strict semantic alignment (ODEs, rules, events, doses, units, constraints) with reproducible, CI‑gated validation suitable for manuscripts and long‑term maintenance.

What’s in here
--------------

- Frozen‑snapshot runner that reproduces SimBiology semantics: `src/offline/frozen_model.py`.
- CLI tools to validate alignment, summarise metrics and plot figures:
  - `python -m scripts.validate_surrogate` — run scenarios, compare against frozen MATLAB references, emit diagnostics.
  - `python -m scripts.summarize_equivalence` — aggregate alignment tables for paper‑ready summaries.
  - `python -m scripts.plot_equivalence` — generate the core figures (boxplots, overlays, convergence).
- Diagnostics: t=0 checks, solver banner logging, per‑scenario event logs (immediate/delayed with trigger time and delay), alignment worst‑point pinpointing.
- An “event suite” snapshot showcasing immediate+delayed events for end‑to‑end testing (`artifacts/matlab_frozen_model/event_suite`).

Repository layout
-----------------

- `src/offline/` — frozen snapshot runner and helpers (`frozen_model.py` is the core).
- `scripts/` — command‑line tools: validation, summarisation, plotting, snapshot/schema utilities.
- `artifacts/matlab_frozen_model/` — frozen snapshots (inputs to the Python runtime).
- `artifacts/validation/` — per‑scenario candidate/reference trajectories and logs.
- `artifacts/analysis/` — aggregated CSV tables for publication.
- `plots/` — generated figures.
- `docs/` — semantics and analysis plan.

What changed vs. MATLAB / original QSP‑IO
------------------------------------------

- Single source of truth: `equations.txt` from `getequations(model)` — Python evaluates the RHS exactly; no S·v reconstruction.
- Rules and constraints
  - Repeated assignments are evaluated in topological order at every RHS and after events.
  - Rate rules overlay on ODE derivatives; algebraic rules solved symbolically where possible, otherwise residuals are asserted.
  - Species flags respected: `BoundaryCondition`, `Constant`, `NonNegative` (projection at sync time).
- Events and doses
  - Triggers support direction; same‑time events execute in SimBiology order; delays are queued and applied at `t_event+ε`.
  - Doses map to events (bolus) or exogenous rates (infusion); t=0 doses precede integration.
- Units and time
  - Honour `configset.json` for time units and tolerances. Compartment volumes and species dimensions are interpreted from snapshot metadata; explicit conversions only where needed.
- Engineering hardening
  - t=0 diagnostic gate (fluxes/ODE RHS/triggers) + trajectory alignment gate (default `max_rel_err ≤ 1e-6`).
  - SHA digests embedded in references for provenance; banner logs solver config and snapshot fingerprints.
  - CI workflow runs tests, alignment, summary and plotting and uploads diagnostics on failure.

Install
-------

Requirements: Python 3.10+, NumPy, SciPy, Pandas, SymPy, Matplotlib, PyTest, Pydantic.

```bash
python -m pip install -U pip
python -m pip install numpy scipy pandas sympy matplotlib pytest pydantic
```

Quick start
-----------

1) Validate scenarios and emit diagnostics

```bash
python -m scripts.validate_surrogate \
  --output artifacts/validation \
  --scenarios all \
  --seed 12345 \
  --emit-diagnostics \
  --max-rel-err 1e-6
```

By default this compares frozen Python semantics (`simulate_frozen_model`) to MATLAB snapshot references for each scenario in the registry.

2) Summarise alignment

```bash
python -m scripts.summarize_equivalence \
  --validation-dir artifacts/validation \
  --output-dir artifacts/analysis
```

Outputs under `artifacts/analysis/`:

- `alignment_metrics.csv` — per‑scenario/per‑observable max relative error, rRMSE, RMSE, R², ΔAUC.
- `alignment_summary.csv` — per‑scenario max rel error and rRMSE quantiles.
- `model_scale.csv` — species/reactions/events/rules/doses counts per snapshot.
- `event_logs.csv` — concatenated per‑scenario event logs with residuals.
- `event_summary.csv` — per‑scenario event timing residual stats.

3) Produce figures

```bash
python -m scripts.plot_equivalence \
  --analysis-dir artifacts/analysis \
  --validation-dir artifacts/validation \
  --plots-dir plots
```

Generated under `plots/`:

- `max_rel_err_boxplot.png` — distribution of max rel error per scenario.
- `event_time_residuals.png` — event timing residuals per scenario.
- `overlay_<scenario>_<observable>.png` — reference vs Python overlays for selected observables.
- `convergence_<scenario>.png` — tolerance sweep (log–log) vs a tight baseline.

Snapshot directory layout
-------------------------

`load_frozen_model` accepts either a snapshot name under the default
`artifacts/matlab_frozen_model/` directory or an explicit path to a directory
containing the frozen SimBiology export. Each snapshot folder must include at
least the following files:

- `configset.json`
- `equations.txt`
- `species.csv`
- `parameters.csv`
- `compartments.csv`
- `rules.csv`
- `reactions.csv`
- `events.csv`
- `doses.csv`
- `stoichiometry.csv`

Optional metadata such as `variants.csv`, `equations_eval_t0_reference.csv`, or
diagnostic artefacts will be picked up automatically when present. Example:

```python
from pathlib import Path
from src.offline import load_frozen_model

snapshot_dir = Path("/data/qspio_snapshots/tumour_cohort_v2")
model = load_frozen_model(snapshot_dir)
```

Minimal example
----------------

```python
from pathlib import Path

from src.offline import load_frozen_model, simulate_frozen_model, DoseEntry

model = load_frozen_model("example1")
custom_dose = DoseEntry(
    index=9001,
    name="one_off",
    dose_type="RepeatDose",
    target="V_C.nivolumab",
    amount=200.0 / 1.436e8,
    amount_units="mole",
    start_time=0.0,
    interval=1e9,
    repeat_count=0,
)
result = simulate_frozen_model(
    "example1",
    days=28.0,
    therapy="anti_pd1",
    custom_doses=[custom_dose],
    sample_interval_hours=12.0,
)
output_path = Path("artifacts/debug/example1_oneoff.csv")
result.save_csv(output_path)
print(result.to_frame(order="contract").tail())
```

Extra output plugin (AUC demo)
------------------------------

```python
import numpy as np

from src.offline.simulation import ExtraOutputs

class AUCPlugin(ExtraOutputs):
    name = "auc_pd1"

    def compute(self, time_days, states, contexts, meta):
        occupancy = np.array([ctx.get("H_PD1_C1", 0.0) for ctx in contexts])
        auc = np.trapezoid(occupancy, time_days)
        return {"pd1_auc": np.full_like(time_days, auc)}

result = simulate_frozen_model(
    "example1",
    days=14.0,
    therapy="anti_pd1",
    extra_outputs=[AUCPlugin()],
)
```

Programmatic use
----------------

```python
from src.offline.frozen_model import simulate_frozen_model

# Run a frozen snapshot (e.g. example1) in Python
res = simulate_frozen_model("example1", days=400.0, therapy="anti_pd1", seed=12345, emit_diagnostics=True)
df = res.to_frame()

# Access event logs
event_log = []
res = simulate_frozen_model("event_suite", days=10.0, therapy="none", event_log=event_log)
print(event_log[:2])
```

Interfaces (Inputs / Outputs)
-----------------------------

Python API
- `simulate_frozen_model(snapshot, *, days, therapy, seed=None, emit_diagnostics=False, run_label=None, event_log=None, rtol_override=None, atol_override=None) -> ScenarioResult`
  - Inputs:
    - `snapshot` (str): name under `artifacts/matlab_frozen_model/` (e.g. `example1`, `event_suite`).
    - `days` (float): simulation horizon in `configset.TimeUnits` (days by default).
    - `therapy` ("none" | "anti_pd1"): treatment mode used by the snapshot.
    - `seed` (int|None): deterministic seed (banner logs include it).
    - `emit_diagnostics` (bool): if true, prints solver banner and enables event logging.
    - `event_log` (list|None): if provided, populated with dict items per firing (see schema below).
    - `rtol_override`/`atol_override` (float|None): temporarily override solver tolerances (for sweeps).
  - Output: `ScenarioResult` with arrays
    - `time_days`, `cancer_cells`, `dead_cells`, `t_cells`, `tumour_volume_l`, `tumour_diameter_cm`, `pd1_occupancy`, `tcell_density_per_ul`.
    - `ScenarioResult.to_frame()` returns a DataFrame with the same columns.
  - Event log entry (when `event_log` is provided):
    - `{ "event_index": int, "time_fire": float, "time_trigger": float, "delay": float, "type": "immediate"|"delayed", "assignments": "Target=Expr; …" }`.

CLI — validation
- `python -m scripts.validate_surrogate [--output PATH] [--scenarios NAMES|all] [--seed INT] [--max-rel-err FLOAT] [--emit-diagnostics]`
  - Inputs:
    - Scenarios come from the registry in `scripts/validate_surrogate.py` (e.g. `example1_control`, `example1_treated`, `example2_treated`, `event_suite`).
    - Frozen snapshot references must exist under `artifacts/validation/*_reference.csv`.
  - Outputs (under `--output`, default `artifacts/validation`):
    - `*_reference.csv`, `*_surrogate.csv` — paired trajectories (reference vs Python candidate).
    - `events_<scenario>_python.csv` — per‑scenario event log (schema below; may be empty but will have headers).
    - `metrics.csv` — per‑observable alignment metrics.
    - `performance.json` — timing snapshot of a reference and surrogate run (keys: `reference_seconds`, `reference_replicates`, `surrogate_seconds`, `surrogate_replicates`).
    - `artefacts.csv` — registry of generated files.
  - Exit: non‑zero if `max_rel_err` gate fails or required files are missing.

CLI — summarisation
- `python -m scripts.summarize_equivalence --validation-dir artifacts/validation --output-dir artifacts/analysis`
  - Inputs: per‑scenario CSVs written by the validator.
  - Outputs (all CSV):
    - `alignment_metrics.csv` — columns: `scenario, observable, max_fractional_error, relative_rmse, rmse, r2, delta_auc`.
    - `alignment_summary.csv` — per scenario: `scenario, max_rel_err, median_rel_rmSE, p95_rel_rmse`.
    - `model_scale.csv` — per snapshot: `snapshot, species, reactions, events, rules, doses`.
    - `event_logs.csv` — concatenated `events_<scenario>_python.csv` plus `time_residual=fire-trigger` and `delay_residual=time_residual-delay`.
    - `event_summary.csv` — `scenario, max_time_residual, median_time_residual, max_delay_residual, median_delay_residual`.

CLI — plotting
- `python -m scripts.plot_equivalence --analysis-dir artifacts/analysis --validation-dir artifacts/validation --plots-dir plots [--scenarios …] [--convergence-scenario NAME] [--seed INT] [--tolerances …]`
  - Inputs: CSV tables from `artifacts/analysis` and trajectories from `artifacts/validation`.
  - Outputs (PNG files under `plots/`):
    - `max_rel_err_boxplot.png`, `event_time_residuals.png`.
    - `overlay_<scenario>_<observable>.png` for each requested scenario.
    - `convergence_<scenario>.png` computed via tolerance sweep against a tight baseline.
  - Fails fast (non‑zero) if required inputs are missing or tables are empty.

CLI — snapshot schema validation
- `python -m scripts.validate_snapshot artifacts/matlab_frozen_model/<snapshot>`
  - Validates file presence and row‑level types using Pydantic; exits non‑zero on the first failure with an informative message (filename and row number).

CSV schemas (quick reference)
- Trajectories: columns `time_days`, `cancer_cells`, `dead_cells`, `t_cells`, `tumour_volume_l`, `tumour_diameter_cm`, `pd1_occupancy`, `tcell_density_per_ul`.
- Event logs (per scenario): `scenario, event_index, time_fire, time_trigger, delay, type, assignments`.
- Metrics: `scenario, observable, rmse, r2, delta_auc, relative_rmse, max_fractional_error`.
- Model scale: `snapshot, species, reactions, events, rules, doses`.
- Event summary: `scenario, max_time_residual, median_time_residual, max_delay_residual, median_delay_residual`.

Event log schema (diagnostics)
------------------------------

When diagnostics are enabled the runtime writes `events_<scenario>_python.csv` with columns:

- `scenario` — scenario id; `event_index` — model order; `time_fire` — execution time; `time_trigger` — first trigger; `delay` — evaluated delay; `type` — `immediate|delayed`; `assignments` — canonical `Target=Expr` string.

Validation results (current snapshots)
--------------------------------------

- Example 1/2 and the event suite align within numerical noise (max rel error ~ 10⁻¹⁵ for tracked observables). The CI gate remains at `max_rel_err ≤ 1e-6` with `RelTol=1e-6, AbsTol=1e-12`.
- Event timing residuals for the event suite reflect the designed immediate (`delay=0`) and delayed (`delay=1`) firings; residual summaries are exported under `artifacts/analysis/event_summary.csv`.

Exporting snapshots from MATLAB
-------------------------------

Use `matlab/scripts/export_matlab_snapshot.m` to freeze a SimBiology project (assumes project scripts akin to `matlab/scripts/example1.m`). The exporter writes `equations.txt`, `configset.json`, `species.csv`, `parameters.csv`, `compartments.csv`, `rules.csv`, `events.csv`, `doses.csv`, `variants.csv`, `reactions.csv`, `stoichiometry.csv` under `artifacts/matlab_frozen_model/<snapshot>/`.

> The MATLAB helpers are not tracked in this Python-only repository. Retrieve them from the original QSPIO distribution (see References below) and keep them locally under the git-ignored `matlab/` directory.

Adding a new scenario (step‑by‑step)
------------------------------------

1) Export or place a snapshot under `artifacts/matlab_frozen_model/<name>/` (see above).
2) Register the scenario in `scripts/validate_surrogate.py` (`SCENARIO_REGISTRY`) with:
   - a scenario `name` (e.g., `my_model_control`),
   - the parameter file tuple (for bookkeeping),
   - `therapy` ("none" or "anti_pd1"),
   - `snapshot` (must match the folder under `artifacts/matlab_frozen_model/`),
   - optional `stop_time` (defaults to 400.0 days).
3) Generate the MATLAB reference CSV once (you can temporarily use the Python reproduction to bootstrap it):
   ```python
   from pathlib import Path
   from src.offline.frozen_model import simulate_frozen_model
   df = simulate_frozen_model("<snapshot>", days=400.0, therapy="none").to_frame()
   Path("artifacts/validation").mkdir(parents=True, exist_ok=True)
   df.to_csv("artifacts/validation/<scenario>_reference.csv", index=False)
   ```
   Replace this file with the true MATLAB trajectory if/when available.
4) Run validation/summarisation/plotting as in the Quick start.

Troubleshooting
---------------

- Missing reference CSV — create `artifacts/validation/<scenario>_reference.csv` (see “Adding a new scenario”).
- Empty event logs — expected for models without event firings; the CSV will still have headers.
- Gate failure (`ALIGNMENT_FAIL …`) — use `--emit-diagnostics` to print solver banners and inspect per‑scenario metrics/logs under `artifacts/validation/`.
- Schema errors — run `python -m scripts.validate_snapshot artifacts/matlab_frozen_model/<snapshot>` and fix the reported row/column.
- Convergence — use `scripts.plot_equivalence --convergence-scenario <name> --tolerances …` to spot ill‑conditioned models.

Integrating QSPIO‑TNBC (planned)
---------------------------------

The TNBC extension can be brought into this Python runtime by freezing its SimBiology project and registering the resulting snapshots as scenarios.

Suggested steps:

1) Export TNBC snapshots (e.g., `tnbc_control`, `tnbc_treated`) using `matlab/scripts/export_matlab_snapshot.m` and place them under `artifacts/matlab_frozen_model/tnbc_*`.
2) Register the scenarios in `SCENARIO_REGISTRY` with appropriate `therapy` and `stop_time`.
3) Generate initial reference CSVs in `artifacts/validation/tnbc_*_reference.csv` (from MATLAB when available, else bootstrap with Python reproduction and replace later).
4) Validate and summarise:
   ```bash
   python -m scripts.validate_surrogate --scenarios tnbc_control tnbc_treated --emit-diagnostics --seed 12345 --max-rel-err 1e-6
   python -m scripts.summarize_equivalence --output-dir artifacts/analysis
   python -m scripts.plot_equivalence --plots-dir plots --scenarios tnbc_control tnbc_treated
   ```
5) If TNBC introduces additional observables, they will automatically appear in the overlays; alignment metrics are computed per observable.

References
----------

- Popel Lab, “QSPIO” GitHub repository — https://github.com/popellab/qspio
- Richard Sové, “QSPIO Version 1.0” MATLAB Central File Exchange entry — https://www.mathworks.com/matlabcentral/fileexchange/77234-qspio-version-1-0
- R. J. Sové et al., “QSP-IO: A Quantitative Systems Pharmacology Toolbox for Mechanistic Multiscale Modeling for Immuno-Oncology Applications.” *CPT: Pharmacometrics & Systems Pharmacology* 9(9):484-497, 2020. doi:10.1002/psp4.12546. (“The code can be downloaded from our GitHub page (www.github.com/popellab/qspio).”)

Exporting snapshots from MATLAB
-------------------------------

Use `matlab/scripts/export_matlab_snapshot.m` to freeze a SimBiology project (assumes project scripts akin to `matlab/scripts/example1.m`). The exporter writes `equations.txt`, `configset.json`, `species.csv`, `parameters.csv`, `compartments.csv`, `rules.csv`, `events.csv`, `doses.csv`, `variants.csv`, `reactions.csv`, `stoichiometry.csv` under `artifacts/matlab_frozen_model/<snapshot>/`.

Continuous Integration
----------------------

The GitHub Actions workflow installs dependencies, runs the test suite, validates alignment with diagnostics, summarises metrics, generates figures, and uploads diagnostics on failure (and analysis assets on success). This makes the semantic alignment claims repeatable across platforms.

Limitations
-----------

- Focused on ODE models; no SSA/PDE support. Algebraic solving uses SymPy and falls back to residual checks when closed forms cannot be found. Snapshot schema must match the exporter (validated by `scripts/validate_snapshot.py`).

Contributing
------------

- Add new snapshots under `artifacts/matlab_frozen_model/<name>` and register scenarios in `scripts/validate_surrogate.py`.
- Keep validation references (`artifacts/validation/*_reference.csv`) up to date after exporter/model changes.
- Extend tests for events/units/solver options as needed; PRs welcome.

License
-------

See repository metadata; third‑party MATLAB models remain under their original terms.
