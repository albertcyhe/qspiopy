# Semantic Equivalence of SimBiology Snapshots in Python

This document sketches the analysis plan and manuscript outline for a paper describing the frozen SimBiology → Python execution pipeline. The aim is to provide repeatable, quantitative evidence that the Python runtime preserves SimBiology semantics and to position the system as a practical, CI-friendly alternative when MATLAB is unavailable.

## 1. Claims

1. **Semantic equivalence** – for the supported subset (ODEs, unit conversions, repeated assignments, rate/algebraic rules, events, doses, constraints) the Python integrator reproduces MATLAB/SimBiology trajectories within solver tolerances.
2. **Engineering reproducibility** – schema validation, t=0 diagnostics, and trajectory alignment act as automated gates in CI and succeed consistently across platforms.
3. **Generalisability** – frozen snapshots generated from heterogeneous SimBiology projects (including third-party SBML imports) satisfy the same equivalence criteria without code changes.

## 2. Experimental Design

### 2.1 Model Suite

| Model | Features | Purpose |
|-------|----------|---------|
| QSP-IO Example 1 | Events, doses, RA topology, unit conversion | Baseline tumour-immune case |
| QSP-IO Example 2 | Multiple concurrent events, multi-compartment | Baseline treated control |
| PK 1-/2-compartment + infusion | Infusion start/stop, unit handling | Dose scheduling validation |
| Feedback loop with repeated assignments | Topological evaluation order | RA stability |
| Algebraic rule steady-state subsystem | Symbolic algebraic solver | Algebraic consistency |
| Multi-event with delays | Delay queue, event ordering | Temporal semantics |
| External SBML-derived model(s) | Third-party validation | Generalisability |

### 2.2 Metrics

For each observable \(k\) sampled at times \(t_i\):

* **Maximum relative error** \( e_{\infty}(k) = \max_i \frac{|y^{\text{py}}_k(t_i) - y^{\text{mat}}_k(t_i)|}{\max(|y^{\text{mat}}_k(t_i)|, \epsilon)} \) with \(\epsilon = 10^{-12}\).
* **Relative RMSE** \( \text{rRMSE}(k) = \sqrt{\frac{1}{N} \sum_i \big( \frac{y^{\text{py}}_k(t_i) - y^{\text{mat}}_k(t_i)}{\max(|y^{\text{mat}}_k(t_i)|, \epsilon)} \big)^2 } \).
* **Event time residuals**: \( \Delta t_e = |t^{\text{py}}_e - t^{\text{mat}}_e| \); report max, median, and distribution.
* **Constraint violations**: number of NonNegative projections and Boundary/Constant breaches (expect 0).
* **t=0 diagnostic parity**: boolean consistency across fluxes, ODE RHS, event triggers.

### 2.3 Acceptance Thresholds

* Alignment gate set at `max_rel_err ≤ 1e-6`, matching solver tolerance (`RelTol=1e-6`, `AbsTol=1e-12`).
* Event time differences ≤ `1e-9 × StopTime` in the configured time units.
* No constraint violations (log assertions otherwise).

### 2.4 Statistical Products

* **Table:** model summary (species, reactions, events, rules counts).
* **Table:** per-model max relative error, rRMSE, event time statistics.
* **Box/violin plots:** distribution of max relative errors across models.
* **Histogram/box plot:** event time residuals.
* **Line overlays:** representative observables (Python vs MATLAB) for key models.

## 3. Ablation & Robustness

* **Tolerance sweep:** vary `RelTol/AbsTol` (1e-3…1e-7) and show log–log convergence of alignment error.
* **Semantic ablations:** disable events, repeated assignments, or unit conversions to quantify resulting errors and demonstrate necessity.
* **Cross-solver/platform:** compare results using distinct solvers (MATLAB: `ode15s`, `ode23t`; Python: BDF, LSODA) and OS targets (macOS, Linux). Report residual differences.

## 4. Performance Characterisation (Optional)

* Wall-clock timings and memory footprint for each model under both runtimes.
* Scaling with respect to model size (number of ODE states/events).
* Discuss trade-offs—goal is reproducibility rather than absolute speed.

## 5. Methodology & Reproducibility

* **Pipeline illustration:** freezing scripts → schema validation (`scripts/validate_snapshot`) → t=0 checks → alignment (`scripts/validate_surrogate`).
* **Artifacts:** embed equations/config SHA digests in reference CSVs; publish scripts, frozen snapshots, CI workflow, and diagnostic outputs.
* **CI evidence:** include workflow excerpt showing schema/t0/alignment gates and artifact uploads.
* **Environment capture:** provide `environment.yml` or Dockerfile; document MATLAB export procedure for new snapshots.

## 6. Limitations

* Scope excludes discrete/SSA, PDE solvers, GUI elements, advanced SimBiology features (variants with rules, SBML events with priorities, etc.).
* Algebraic solution relies on SymPy symbolic solver; note potential limitations for stiff or multi-root systems.
* Future work: SBML import automation, stochastic solvers, support for additional SimBiology constructs.

## 7. Manuscript Outline

1. **Introduction** – motivation for MATLAB-free reproducibility in immuno-oncology QSP; overview of SimBiology semantics.
2. **Methods** – export pipeline, schema checks, Python runtime architecture, diagnostics, CI integration.
3. **Results** – semantic equivalence metrics, ablations, convergence, cross-platform results.
4. **Discussion** – implications for model sharing, limitations, extensions to other domains.
5. **Supplementary** – full metric tables, CI logs, schema definitions, automated scripts.

## 8. Next Steps

1. Freeze additional benchmark models and generate MATLAB reference trajectories.
2. Automate metric computation (`scripts/validate_surrogate --emit-diagnostics`) and aggregate into publishable tables/plots (`python -m scripts.summarize_equivalence`).
3. Extend testing to cover infusion endpoints, delayed events, and algebraic edge cases; record constraint violation logs.
4. Set up CI artifact uploads (alignment summaries, diagnostics) for fail-fast reproducibility evidence.
