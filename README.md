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
