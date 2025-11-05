# TODO – SimBiology semantic alignment plan

> Goal: replay SimBiology execution semantics in a pure‑Python environment (no MATLAB required) and match MATLAB reference trajectories within tight thresholds.

---

## 0. Principles
- `getequations` is the single source of truth: evaluate `Fluxes:`/`ODEs:` RHS directly; do not reconstruct via `S·v`.
- No global unit rewriting: interpret units/dimensions from `configset.json` and object metadata; convert only where semantics require.
- Preserve evaluation order: repeated assignments (each step and after events), concurrent events by model order, and t=0 initialisation order (initial values/variants/initial assignments/doses).

---

## 1. Priorities (P0 → P2)

### P0 (must‑have blockers)
- [x] Honor `configset.json` for **TimeUnits/StopTime/SolverOptions/CompileOptions (UnitConversion/DimensionalAnalysis)**.
- [x] Units/dimensions inference: derive `interpreted_dimension` from `species.csv.units` and `CompileOptions.DefaultSpeciesDimension` (amount / concentration / amount/volume, etc.).
- [x] Volume normalisation (e.g., `V_T`): if repeated assignments accumulate µm³, convert to litres per UnitConversion (e.g., `V_T_L = V_T_um3 * 1e-15`).
- [x] ODE direct evaluation: in `f(t, y)` evaluate repeated assignments → compute `Fluxes:` (cacheable) → evaluate `ODEs:` line by line.
- [x] t=0 quick check: compare `Fluxes:`/`ODEs:` RHS against frozen MATLAB; event trigger booleans must match.

### P1 (core semantics)
- [x] Repeated assignments: build dependency graph from expression AST, evaluate in topological order; re‑evaluate after events.
- [x] Event scheduling: execute all same‑time triggers in model order; re‑evaluate repeated assignments; restart integration at `t_event + ε`.
- [x] Doses: map bolus to immediate events; infusion to exogenous rates across the interval; apply t=0 doses before integration.
- [x] Species flags: enforce `BoundaryCondition`/`Constant`/`NonNegative` (clamp nonnegative after sync; reactions do not change boundary species).

### P2 (polish)
- [x] Rate/Algebraic rules: inject rate rules into derivatives; solve algebraic rules via symbolic solution or residual assertions.
- [x] Diagnostics/regression: export `equations_eval_t0.csv` (variable, python_value, reference_value, rel_err); enforce trajectory RMSE thresholds.
- [ ] Optional oracle: SBML + RoadRunner as a difference magnifier (treat `getequations` as ground truth).

---

## 2. Export contract

### 2.1 `configset.json`
- `TimeUnits`, `StopTime`
- `SolverType`, `SolverOptions`（`RelTol`, `AbsTol`, `MaxStep`）
- `CompileOptions`（`UnitConversion`, `DimensionalAnalysis`, `DefaultSpeciesDimension`）

### 2.2 `species.csv` (augmented columns)
- `name`, `compartment`, `initial_value`, `initial_units`
- `units` (string)
- `boundary_condition` (bool), `constant` (bool), `nonnegative` (bool, optional)
- `interpreted_dimension` (derived: amount / concentration / amount/volume / other)

### 2.3 Others
- `rules.csv`: `type(initial|repeated|rate|algebraic)`, `target`, `expr`
- `events.csv`: `index_in_model`, `trigger`, `delay`, `assignments[]`
- `doses.csv`: `type(bolus|infusion)`, `target`, `amount/units`, `time/interval/repeat`, `rate/units`
- `variants.csv`: variant name, assignment list, active flag
- `equations.txt`: verbatim `getequations(model)` output

---

## 3. Python runtime semantics (frozen_model.py)

### 3.1 ODE evaluation order
1) Parse config and construct the symbol table (t, y, parameters, compartments).  
2) Evaluate repeated assignments to a fixed point (once per step).  
3) Compute and cache **Fluxes**.  
4) Evaluate `ODEs:` line by line to obtain `dy/dt`.  
5) Apply **Rate/Algebraic rules** (when present).  

### 3.2 Events and doses
- Trigger via `solve_ivp` events (or custom root‑finding); execute concurrent events at the same time in `index_in_model` order.  
- Re‑evaluate repeated assignments after applying event assignments; continue integration from `y(t_event+)` with `t ← t_event + ε`.  
- Doses: bolus → instantaneous event; infusion → exogenous rate over the interval (or temporary state/rule).  
- Apply t=0 doses and initial assignments following SimBiology’s initialization order.

### 3.3 Species flags and constraints
- `BoundaryCondition=true`: reactions do not change this species; a rate rule (if present) still governs it.  
- `Constant=true`: disallow any reaction/rule updates during the simulation.  
- `NonNegative=true`: project/clamp to 0 after steps or via event‑style truncation.

---

## 4. Units and dimensions
- Use `configset.TimeUnits` as the time base; do not hand‑rewrite all parameters to per‑day/per‑minute.  
- For volumes: perform conversions only where semantically required (e.g., `V_T_um3 → L` with factor `1e-15`).  
- Concentration vs amount: if `units` are missing, use `DefaultSpeciesDimension`; otherwise infer from unit dimensions (amount, amount/volume, …).  
- Presence of `1/compartment` in `ODEs:` is a strong signal for concentration states; do not rewrite RHS based on this alone.

---

## 5. t=0 quick check (must pass before updating references)
- [x] Compute `Fluxes:` item by item and compare against the frozen reference (relative error ≤ 1e‑9 or configured threshold).  
- [x] Compute `ODEs:` RHS item by item and compare.  
- [x] Evaluate all `events.trigger` booleans and require consistency.  
- [x] Export `equations_eval_t0.csv`: `name, python_value, reference_value, rel_err`.

---

## 6. Definition of Done
- [ ] Example 1–N: RMSE and maximum relative error for key observables are within thresholds (derived from model scale and solver tolerances).  
- [ ] All event timestamps align within tolerance (± event tolerance).  
- [ ] No violations of NonNegative/Boundary/Constant constraints (add assertions and logs).  
- [ ] Unit test coverage: expression parsing, dependency topology, concurrent event ordering, dose application, t=0 quick check.  

---

## 7. Logging and diagnostics (suggested)
- `--trace-step`: print each triggered event and the repeated‑assignment targets/values updated per step.  
- `--trace-units`: print unit inference and key conversion factors (e.g., µm³→L).  
- `--dump-t0`: emit Flux/ODE/event assertions at t=0.

---

## 8. Out of scope (for this phase)
- Global unit system auto‑conversion (avoid double scaling on top of `getequations`).  
- Advanced numerical stability optimisations for complex algebraic systems (focus on correctness first).

---

## 9. Final sprint (next steps)
- [x] CI gating: `python -m scripts.validate_surrogate --replicates 1 --max-rel-err 1e-6`, fail fast on t=0 diagnostics.
- [ ] Reference trajectory versioning: record `equations.txt` / `configset.json` SHAs; require explicit bump to replace.
- [ ] Events/doses regression: add tests for concurrency, delay ε stepping, trigger direction tolerance.
- [ ] Algebraic stress tests: multi‑root/nonlinear cases; validate solver fallback and logging.
- [ ] Constraint monitoring: count and log NonNegative/Boundary/Constant violations.
- [ ] Snapshot schema validation: Pydantic/JSON schema at import with version warnings.
- [ ] API/CLI hardening: unify `simulate_frozen_model` signature; CLI `--stop-on-t0-fail`.
- [ ] Cross‑platform validation: Mac/Linux/Windows runs for error thresholds and event timestamps.
