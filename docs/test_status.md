# QSP‑IO Python vs MATLAB Alignment — Test Log

## Executed

- **Example1 / Example2 frozen snapshot regressions**  
  - Runner: `python -m scripts.validate_surrogate --scenarios example1_control example1_treated example2_treated --emit-diagnostics`  
  - Artefacts: `artifacts/validation/` (reference & surrogate CSVs, metrics, solver logs).  
  - Result: ✅ pass — max relative error ≈7e‑15, confirming the frozen snapshot pipeline still mirrors the MATLAB exports.

- **A1 reference regeneration (MATLAB replay)**  
  - Runner: `python -m scripts.run_alignment_suite --scenarios A1 --output artifacts/validation`  
  - Artefacts:  
    - `A1_reference.csv`, `A1_reference_dose_audit.csv`, `A1_events_reference.csv` under `artifacts/validation/`.  
    - Matching surrogate outputs (`A1_surrogate.csv`, header manifest, event log, dose audit).  
  - Result: MATLAB batch run (via `/Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab`) completed and produced the refreshed reference set; post-run metrics step still needs polishing but all CSVs required by the gate now exist.

- **A1 numeric gate (`validate_surrogate`)**  
  - Runner: `python -m scripts.validate_surrogate --scenarios A1 --dump-t0 --numeric-gates --ic-mode snapshot`  
  - Artefacts: `artifacts/validation/A1_*.csv` plus t₀ diagnostic tables inside the same directory.  
  - Result: ❌ fail — PD‑1 synapse转换已生效（`pd1_bridge_block` 将 mol/L 转成 molecule/µm²，默认深度 1.15×10⁻⁵ µm），`pd1_occupancy` 振幅恢复到 MATLAB 量级 (~0.126)，但波形依旧偏平；`tumour_volume_l`/`tcell_density_per_ul` 仍滞于 1e-11 L / 1e-6 cells/µL 量级。该剩余差异被纳入 M4 backlog（动态体积/多克隆），不再阻塞 M3 的交付。
- **Runtime plumbing + 状态事件去抖（M3 部分）**  
  - Changes: `simulate_frozen_model` 默认 `ic_mode="target_volume"`（0.5 cm、150 d、`reset_policy="cancer_only"`；不支持的快照自动回退至 snapshot 初值），并在 `segment_integrator` 中实现 state-trigger 事件的 armed/disarmed 去抖（不应期 + 回到 false 再触发 + ε‑bump）。`scripts/validate_surrogate` CLI 暴露 `--ic-mode/--ic-target-diam-cm/--ic-max-days/--ic-max-wall-seconds/...` 以保持 runner/CI 一致。  
  - Tests: `pytest tests/test_alias_injection.py`、`pytest tests/test_initial_conditions.py -m slow`、`pytest tests/test_events.py -k event_suite`、`pytest tests/test_state_trigger_hysteresis.py`。  
  - Status: ✅ 防抖语义已在轻量诊断场景通过；尚未回到 example2/A1 snapshot 模式做整体验证，因而 runner 仍以 target-volume IC 为默认，待 snapshot 路径验证后再切回。
- **Warm-start kick + snapshot example2 复测（M3 补丁）**  
  - Changes: 在 `segment_integrator` 中增加 deterministic kick `_kick_off_t0()` 与新的 `warm_start_quarantine()`（Heun 推离 + Radau/BDF 补齐 + 墙钟/重试护栏 + `reconcile(vec, time)`），并在 `pd1_bridge_block` 中把 nivolumab PK 自动转成 synapse 面密度（Avogadro × depth），确保 `ic_mode='snapshot'` 不再卡在 t≈0 且 PD‑1 模块具备 2D↔3D 语义。  
  - Tests: `pytest tests/test_state_trigger_hysteresis.py`、`pytest tests/test_events.py -k event_suite`、`pytest tests/test_alias_injection.py`、`python -m scripts.diagnose_t0 example2`。  
  - Result: ✅ `python - <<'PY' ... simulate_frozen_model('example2', ..., ic_mode='snapshot')` 现已能稳定返回，`aPD1_concentration_molar` / `aPD1_surface_molecules_per_um2` 也会随 context 输出。后续的 tumour/occupancy/tcell 波形对齐已转交 M4。

## In Progress

- **A1 semantic debugging** — reconcile tumour volume / PD‑1 occupancy / T-cell density between MATLAB replay and Python surrogate (targets: rel_L2 ≤ 1e‑3, maxRE ≤ 5e‑3). Focus now shifts to biology semantics (tumour module, receptor occupancy rules, repeated assignments) since unit conversions and dose audits match.
- **Scenario orchestration** — registry + CLI plumbing complete; M3 switches are available via CLI, but the target-volume IC still needs a calibrated tumour-growth preset (currently stuck at ~0.012 cm). Once this is tuned, re-run `python -m scripts.validate_surrogate --scenarios A1 --ic-mode target_volume ...` to measure the improvement before unlocking the full A-series.

## Planned Runs

- **Dose matrix (A2–A6)** — once A1 aligns, extend the runner to the remaining monotherapy regimens and regenerate metrics.
- **Combination therapy (B1–B3)** — introduce PD‑1 + CTLA‑4 dosing after confirming combination species/parameters exist in the frozen snapshot.
- **Population grid (C)** — implement cohort sampler with shared RNG seeds and summarise ORR/DCR deltas vs MATLAB.
- **Sensitivity sweep (D)** — ±20 % one‑at‑a‑time perturbations with Spearman ranking comparison.
- **Steady-state/guardrail checks (E/F)** — zero‑dose stability, high‑dose robustness and pathway knock‑out diagnostics.

These queued scenarios remain pending until the baseline MATLAB ↔ Python parity (≤1e‑3 relative L2 on tracked states) is restored.
