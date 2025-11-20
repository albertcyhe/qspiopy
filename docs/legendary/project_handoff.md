# QSP‑IO Snapshot Alignment — Handoff Summary (2025‑11‑13)

This document captures the current state of the SimBiology ↔ Python alignment effort so a new contributor can take over without digging through the whole history.

---

## 1. Context & Goal

- **What we ship**: a Python surrogate that integrates frozen SimBiology snapshots (`artifacts/matlab_frozen_model/<snapshot>`) and reproduces MATLAB reference trajectories (`artifacts/validation/<scenario>_reference.csv`).
- **Why it hurts**: the base solver (M1–M3) is stable and unit‑consistent, but A‑series scenarios still fail the numeric gates because PD‑1 occupancy, tumour volume, and t‑cell density stay “flat” in Python while MATLAB shows slow dynamics.
- **Strategy**: keep the snapshot pipeline untouched (“white box”), add a thin “alignment layer” (runtime modules) only where MATLAB semantics are missing, and gradually retire the grey‑box logic as soon as we can export the true equations.

---

## 2. What’s done (kept and tested)

| Area | Status | Notes |
| --- | --- | --- |
| **Integrator & events (M1)** | ✅ | Segmented integration with ε‑bump, trigger specs, pre/post records (`src/offline/segment_integrator.py`, `simulate_frozen_model`). |
| **Units & parameters (M2)** | ✅ | `src/offline/units.py` normalises everything to day/L/mol/M; dose application goes through `normalise_dose_to_species`. |
| **Snapshot runtime (M3)** | ✅ | `simulate_frozen_model` handles snapshot/target‑volume ICs, event de-bounce, warm starts, CLI toggles; module system in `src/offline/modules/switches.py`. |
| **Instrumentation** | ✅ | `--dump-flat-debug` now captures module internals (PD‑1 filter, geometry stats). New script `scripts/dump_snapshot_semantics.py` mirrors MATLAB’s `show_alignment_drivers.m`. |
| **Grey-box alignment layer** | 🚧 | `alignment_driver_block` introduces a minimal PK/PD‑1/volume ODE, but parameters are placeholder values and waveforms still miss MATLAB shapes. |

Artifacts to know:
- `artifacts/show_alignment_example1.txt`: MATLAB report via `matlab/scripts/show_alignment_drivers.m`.
- `artifacts/show_snapshot_example1_filtered.txt`: matching snapshot dump (filtered for PD‑1/geometry keywords).

---

## 3. Known gaps (must fix before gates turn green)

1. **Exporter parity**: MATLAB exporter currently emits only a subset of reactions/rules in its summary. We need to ensure `reactions.csv`/`equations.txt` contain **every** reaction (including PD‑1 transport/binding) so `FrozenModel.rhs()` stays authoritative.
2. **Runtime overrides**: historical modules (`pd1_bridge_block`, `tumour_geometry_*`) disable the very rules we just exported. For A‑series we should run with **only** `alignment_driver_block` until white-box semantics are restored.
3. **Alignment driver calibration**: the new ODEs (PD‑1 binding + logistic tumour volume) are wired up but unfit—occupancy jumps to 0.13 instantly, volume reacts on day 1, and `tcell_density_per_ul` is off by ~1e6×. We need to fit the explicit parameters (`pd1_occ_*`, `geom_*`) against MATLAB references.
4. **Snapshot parameter coverage**: MATLAB `parameters.csv` now includes the alignment knobs (e.g. `pd1_occ_kon_scale`, `geom_growth_per_day`). Once the white-box exporter is fixed, these should move back into SimBiology so Python no longer hardcodes defaults.

---

## 4. Immediate next steps for a new contributor
直接说结论：**你现在这条路（先把 snapshot A1–A6/B 跑通，再考虑 white‑box 重写）是合理的，而且比一上来大改 ODE 风险小很多。**

下面我给你一版可以几乎原样贴进 handoff 文档的 `## 4. Immediate next steps for a new contributor`，按照“接手第一周能干什么”来拆分。

---

## 4. Immediate next steps for a new contributor

> 目标：**在不改 ODE / 规则的前提下**，用现有冻结 snapshot + alignment driver，把 A1 的数值门先拉到“可接受”，为后续扩展到 A2–A6/B 和 white‑box 重写打基础。

### 4.1 熟悉工程地形 & 复现当前基线（Day 0–0.5）

1. 本地把测试跑一遍（确认环境 OK）：

   ```bash
   pytest tests/test_units.py tests/test_module_blocks.py
   pytest tests/test_simulation_sampling.py
   ```

   这些用到 `FrozenModel` 的加载、单位换算和 repeated assignment 流程，是后面所有工作的基础。

2. 跑一遍 **当前的 A1 snapshot 基线**：

   ```bash
   python -m scripts.validate_surrogate \
     --scenarios A1 \
     --ic-mode snapshot \
     --module-block alignment_driver_block \
     --emit-diagnostics \
     --dump-flat-debug 5 \
     --numeric-gates
   ```

   * 把输出目录 `artifacts/validation/` 留好，尤其是：

     * `A1_reference.csv`（MATLAB 轨迹）
     * `A1_surrogate.csv`（Python 轨迹）
     * `A1_*.metrics.csv`（数值门）
     * `A1_flat_debug_*.csv`（alignment driver 内部信号）

3. 快速浏览一遍代码骨架（只看，不改）：

   * `src/offline/snapshot.py::FrozenModel`：了解 snapshot 如何变成 Python 里的 ODE + rules + events。
   * `src/offline/units.py`：单位统一到 day/L/mol 的逻辑。
   * `src/offline/simulation.py::simulate_frozen_model`：把 FrozenModel + doses + events + module_blocks 串起来的地方。
   * `src/offline/modules/switches.py::alignment_driver_block`：灰箱 PD‑1 + tumour follower 逻辑。
   * `scripts/validate_surrogate.py`：A1 场景 registry、CLI 参数、数值门。

**产出**：一份简短的“接手笔记”，记录当前 A1 的三个关键指标（tumour_volume_l / tcell_density_per_ul / pd1_occupancy 的 rel_L2/maxRE）和对应命令。

---

### 4.2 对齐构架：搞清楚灰箱 driver 的“管线”（Day 0.5–1）

目标：新同事要清楚 **信号从哪来、到哪去**，否则后面调参只是在瞎拧。

1. 画一张小图（可以放进 docs/）：

   * 纵向：`dose → alignment_driver_block PK surrogate → pd1_occupancy surrogate → ScenarioResult.observables`
   * 横向：`FrozenModel.doses` → `simulate_frozen_model(... _active_dose_schedule ...)` → `alignment_driver_block` 写 `context['pd1_alignment_*']` & `context['tumour_volume_l']` → `ScenarioResult` 采样这些键。

2. 用 `--dump-flat-debug` 验证这条管线确实在动：

   * 在 `A1_flat_debug_*.csv` 中检查：

     * `pd1_alignment_pk_state`, `pd1_alignment_input`, `pd1_occupancy`
     * `tumour_volume_l`, `tcell_density_per_ul`
   * 画一张简单的 overlay（MATLAB vs Python）：

     * x 轴：time_days
     * y 轴：上述几个 observable + alignment 内部状态

**产出**：一张小图（或者几张 PNG），能让别人一眼看到“Python 现在的波形长什么样、跟 MATLAB 差在哪”。

---

### 4.3 聚焦 PD‑1 波形：先让占有率“像样”（Day 1–2）

现在 alignment driver 已经产生非平坦的信号了，但 `pd1_occupancy` 的形状差很多。下一步只动 **PD‑1 部分，不动几何**。

1. 编一个最小的 **PD‑1 sandbox 脚本**（建议放在 `scripts/dev_pd1_sandbox.py`）：

   * 只跑 A1，强制开启 `alignment_driver_block`，并：

     * 从 `ScenarioResult.to_frame()` 里提取 `pd1_occupancy` & `pd1_alignment_*` 列。
     * 从 `A1_reference.csv` 读出 MATLAB 的 `pd1_occupancy`。
     * 画 overlay（保存到 `artifacts/dev/`）。

2. 把可调参数集中起来：

   * 例如：`pd1_pk_dose_scale`, `pd1_pk_half_max`, `pd1_occ_tau1`, `pd1_occ_tau2`, `pd1_occ_delay`, 以及任何 hill 相关参数。
   * 确认它们都来自：

     * `parameters/example1_parameters.json`
     * `artifacts/matlab_frozen_model/example1/parameters.csv`
       而不是 CLI 临时 override，确保未来快照导出后仍然生效。

3. 用 `scripts.fit_observables.py` 做 **只针对 PD‑1 的单目标拟合**：

   * 先只拟合 `pd1_occupancy`，冻结几何相关权重：

     ```bash
     python -m scripts.fit_observables \
       --scenario A1 \
       --target pd1_occupancy \
       --module-block alignment_driver_block \
       --max-evals 50 \
       --output artifacts/validation/A1_pd1_fit.json
     ```
   * 拟合完成后：

     ```bash
     python -m scripts.validate_surrogate \
       --scenarios A1 \
       --ic-mode snapshot \
       --module-block alignment_driver_block \
       --param-override @artifacts/validation/A1_pd1_fit.json \
       --dump-flat-debug 5
     ```

4. 验收标准（先给一个“软目标”，便于新同事有方向）：

   * `pd1_occupancy`：

     * `rel_L2 < 0.1`
     * `maxRE < 0.3`
   * 如果达不到，就回 PD‑1 sandbox 里看：

     * 是不是 PK surrogate 太快（τ 太小）？
     * 是否需要增加一个 **显式 delay / 二阶 filter**（这时候可以改 alignment_driver_block 的 ODE 形式，而不是再去硬调 τ）。

**产出**：

* `artifacts/validation/A1_pd1_fit.json`（拟合出来的一组 PD‑1 参数）
* 一段简短说明记录目前 `rel_L2/maxRE` 数值 & 使用的配置（方便以后回溯）。

---

### 4.4 再看几何：让体积/密度先“对趋势”（Day 2–3）

在 PD‑1 有个像样的波形后，再轮到 tumour volume & T‑cell density。

1. 同样做一个 **几何 sandbox**：

   * 脚本 `scripts/dev_geometry_sandbox.py`：

     * 读 `A1_reference.csv` 的 `tumour_volume_l` 和 `tcell_density_per_ul`。
     * 读 Python 的 `ScenarioResult` 同名列（确保 alignment driver 写入的是这些键）。
     * 画 overlay。

2. 检查 driver 输入是不是真的在动：

   * 看 `A1_flat_debug_*.csv` 中：

     * `geom_driver_volume_raw`（比如 snapshot 的 V_T 或 cell 总量）
     * `geom_driver_volume_smoothed`
     * `geom_driver_tcell_raw` / `geom_driver_tcell_smoothed`
   * 如果连“driver_raw”都几乎水平，那就说明 ODE 本体没有在变化，这时候：

     * 你可以让 alignment driver 用 **参考 CSV** 做 follower（例如：一阶滞后跟随 reference volume），或者
     * 用一个简单的 “growth + clearance” surrogate：`dV/dt = k_grow * V - k_clear * V`，参数完全独立于 snapshot。

3. 用 `scripts.fit_observables.py` 做 **仅几何相关的拟合**：

   * 在已冻结的 PD‑1 参数基础上，拟合几何部分：

     ```bash
     python -m scripts.fit_observables \
       --scenario A1 \
       --target tumour_volume_l tcell_density_per_ul \
       --module-block alignment_driver_block \
       --fix-params pd1_* \
       --max-evals 80 \
       --output artifacts/validation/A1_geom_fit.json
     ```
   * 应用并复测：

     ```bash
     python -m scripts.validate_surrogate \
       --scenarios A1 \
       --ic-mode snapshot \
       --module-block alignment_driver_block \
       --param-override @artifacts/validation/A1_pd1_fit.json \
       --param-override @artifacts/validation/A1_geom_fit.json \
       --numeric-gates \
       --dump-flat-debug 5
     ```

4. 验收目标（Phase 1）：

   * `tumour_volume_l`：`rel_L2 < 0.2`，`maxRE < 0.3`
   * `tcell_density_per_ul`：先争取 `rel_L2 < 0.5`，几何/密度本身就更敏感，可以稍微宽松一点。

**产出**：

* `A1_geom_fit.json`
* 几张体积/密度 overlay 图 + 一份“误差表”。

---

### 4.5 把拟合结果固化到 snapshot / 参数源（Day 3–4）

现在灰箱 driver 基本能跑出“像样”的 A1 波形了，要把这些参数固化，避免全靠 CLI override。

1. 把 PD‑1 + 几何的最终参数整理成一个 **统一 JSON**：

   * `artifacts/validation/A1_alignment_final.json`
   * 结构可以是 `{ "param_name": value, ... }`，方便直接喂给 `--param-override @...`。

2. 更新 canonical 参数源：

   * 修改 `parameters/example1_parameters.json`：把对齐好的值写回去（注意 units 字段不要乱动）。
   * 重新导出 snapshot（如果你已经有 MATLAB pipeline）或直接编辑：

     * `artifacts/matlab_frozen_model/example1/parameters.csv` 相应行。

3. 跑一次 snapshot 校验：

   ```bash
   python -m scripts.validate_snapshot artifacts/matlab_frozen_model/example1
   ```

   确认 schema 仍然过。

4. 再跑一遍 **无 override** 的 A1 验证：

   ```bash
   python -m scripts.validate_surrogate \
     --scenarios A1 \
     --ic-mode snapshot \
     --module-block alignment_driver_block \
     --numeric-gates \
     --dump-flat-debug 5
   ```

   * 检查数值门是否稳定（哪怕还没完全绿，也要记录当前数值）。

**产出**：

* 更新过的 `example1_parameters.json` / `parameters.csv`
* 一份说明：当前 A1 在“零 override”配置下的数值门结果。

---

### 4.6 扩展到 A2–A6/B & 为 white‑box 铺路（Day 4–5+）

当 A1 的 snapshot 路径基本“可用”后，新同事可以开始做两件事：

1. **把 alignment driver 横向推广到其它 scenario**：

   * 在 `scripts/validate_surrogate.py` 的 `SCENARIO_REGISTRY` 里，把 A2–A6/B 也挂上 `alignment_driver_block`。
   * 先只看趋势（不必立刻拟合所有场景），记录每个场景的 PD‑1 / 体积 / 密度波形和门的误差。

2. **为将来的 white‑box 重写记笔记**：

   * 在一个新文档（比如 `docs/whitebox_transition_plan.md`）里写下：

     * 目前 alignment driver 的 ODE 形式（PK + occupancy + geometry）；
     * 哪些状态将来希望迁入真正的 FrozenModel.odes / reactions；
     * 将来想加的“代谢/微环境/贝叶斯优化 hook”应该挂在哪里（例如：在 `FrozenModel.evaluate_reactions` 中插入 meta‑bolic flux，或者在 alignment driver 基础上增加 control inputs）。



---

## 5. Repo map for handoff

| Path | Purpose |
| --- | --- |
| `src/offline/` | Snapshot loader (`snapshot.py`), solver (`simulation.py`), runtime modules (`modules/switches.py`). |
| `scripts/validate_surrogate.py` | CLI entrypoint for scenarios; registers A/B series via `scripts/scenario_registry.py`. |
| `scripts/fit_observables.py` | Grey/white-box parameter fitter (least-squares). |
| `scripts/dump_snapshot_semantics.py` | Snapshot-side dump for PD‑1/geometry semantics. |
| `docs/new_alignment_plan.md` | Living summary of milestones/status (now trimmed; see below). |
| `docs/m4_geometry_status.md` | Historical notes on geometry calibration; keep for reference. |
| `artifacts/matlab_frozen_model/*` | Frozen snapshots (input). |
| `artifacts/validation/*` | MATLAB reference CSVs, Python surrogate outputs, diagnostics. |

---

## 6. Contact checklist

When syncing with MATLAB-side owners, align on:
- Exporter coverage (full reaction/rule dump).  
- Official PD‑1 geometry constants (`geometry_*`, `pd1_occ_*`).  
- Whether PD‑1 occupancy requires additional MATLAB-only modules (if so, expose their math so we can reproduce it in Python).
- 你可以直接调用 /Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab 来调用本地matlab环境，完成一些任务。


Once these pieces are in place, the remaining work is mainly calibration + regression tests. Good luck!
