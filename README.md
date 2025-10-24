QSPIO - A MATLAB SimBiology-Based Platform for QSP Modelling in Immune Oncology

## Cloning the Repository

Click "Clone or Download" button on the main page.

**or**

Using the command line:

`git clone https://github.com/popellab/qspio.git`

---

## Installing the Repository

Add qspio directory with subdirectories to your MATLAB path.

---

## Offline Python CLI for Contract Tests

The repository ships with a lightweight Python surrogate for the MATLAB
SimBiology workflow.  The CLI reads a tidy `EPOCH_FEATURES.csv`, validates the
content against the offline schema, runs the deterministic surrogate, and
emits a reproducible `QSP_OUT.csv` artefact.  This allows contract tests to run
without MATLAB while retaining the original gating and metadata guarantees.

The CLI lives in `scripts/qsp_io_cli.py` and can be invoked directly via
`python -m`:

```bash
python -m scripts.qsp_io_cli \
  --input artifacts/EPOCH_FEATURES.csv \
  --output artifacts/QSP_OUT.csv \
  --seed 42
```

Key behaviours:

* The input CSV must provide the columns `mechanism_id`, `feature`, `time_h`,
  `value`, and `units`.  Missing columns, empty tables, inconsistent mechanism
  assignments, or unit mismatches terminate the run with a non-zero exit code.
* Mechanism identifiers are validated against
  `parameters/MECHANISM_REGISTRY.csv` when present.  When the registry is not
  available the tool falls back to the identifiers described in the offline
  schema.
* The tool records the SHA-256 digest of the input table, the deterministic
  seed, and the solver's `source_version` metadata alongside the simulated
  trajectories.  These fields make it straightforward to diff artefacts during
  contract testing.
* Use `--dry-run` to exercise the full validation/simulation path without
  writing the `QSP_OUT.csv` artefact.  This is helpful for CI jobs that only
  need to confirm the CSV contract.

The default locations for the artefacts (`artifacts/EPOCH_FEATURES.csv` and
`artifacts/QSP_OUT.csv`) mirror the paths consumed by downstream contract
tests, but both can be overridden via the CLI arguments to match bespoke
pipelines.

## Tutorial reproduction helper

The MATLAB tutorial referenced in the publication ships with rich parameter
catalogues (`parameters/example1_parameters.json`, `parameters/example2_parameters.json`).
To explore those scenarios – or external catalogues such as the TNBC project –
without launching MATLAB you can rely on the reduced-order Python surrogate::

    python -m scripts.reproduce_tutorial \
        --parameters parameters/example1_parameters.json \
        --therapy anti_pd1 \
        --output artifacts/example1_python.csv

The CLI will resolve all derived parameters, run the simplified tumour/T-cell
model, and emit a CSV containing tumour volume, diameter, and cell counts over
time.  Multiple `--parameters` flags can be provided to layer project-specific
overrides on top of the base catalogue, mirroring the MATLAB workflow.

## Consistency and performance validation

To contrast the lightweight surrogate against the mechanistic reference
equations, run:

```bash
python -m scripts.validate_surrogate \
  --output artifacts/validation \
  --seed 42 \
  --emit-diagnostics \
  --scenarios all
```

The command generates paired CSV trajectories for the Example 1/2 scenarios by
comparing the frozen-core Python runtime (`simulate_frozen_model`) against the
MATLAB snapshot references, computes RMSE/R²/ΔAUC/TGI statistics, records a
wall-clock comparison, and emits solver/banner diagnostics with alignment
summaries. The resulting
artefacts power the quantitative assessment documented in
`docs/semantic_equivalence_report.md`.

To consolidate the per-scenario diagnostics into publishable tables run:

```bash
python -m scripts.summarize_equivalence \
  --validation-dir artifacts/validation \
  --output-dir artifacts/analysis
```

The generated CSV files (`alignment_metrics.csv`, `alignment_summary.csv`,
`model_scale.csv`, `event_logs.csv`) provide the backbone for the tables outlined in the report.

To render the core figures described in the paper plan:

```bash
python -m scripts.plot_equivalence \
  --analysis-dir artifacts/analysis \
  --validation-dir artifacts/validation \
  --plots-dir plots
```

This command produces the max-relative-error boxplot, event residual chart,
representative overlay plots, and tolerance sweep curves under `plots/`.
Add `--use-tutorial-surrogate` if you need to compare against the reduced
contract-testing surrogate instead of the frozen snapshot.
