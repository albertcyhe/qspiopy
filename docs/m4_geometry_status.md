# M4 T cell & 动态几何白盒化路线图（2025‑11‑20）

PD‑1 子系统已经照 stiff solver 模板完成白盒迁移（probe → diff → clean export）。接下来我们要用同一套流程，把 T cell 激活/增殖反应以及 tumour geometry follower 也搬进白盒世界，最终让 `alignment_driver_block` 只负责 orchestrate，而非承载任何灰盒逻辑。

---

## 1. 全局视角

| 领域 | 当前状态 | 下一步 |
| --- | --- | --- |
| PD‑1 synapse | ✅ 已和 SimBiology ODE 对齐 (`dev_pd1_training_probe.m` + `dev_pd1_probe_diff.py`)，可导出干净训练集 (`export_pd1_clean_training.py`)。 | 作为模板指导 T cell / 几何白盒化。 |
| T cell 驱动 | ⚠️ 基础设施已到位：`dev_tcell_probe.m` + `dev_tcell_probe_diff.py` 可获取/对比 MATLAB 轨迹，`tcell_whitebox.py` 已接入 stiff solver 并加载 snapshot 参数，但首次 diff RMSE（density≈8e10、Tumor_T1≈1e8）仍远超目标。 | 以 probe diff 为 tight loop 调整 Reaction 14–27 RHS 和参数，压低 RMSE，再接入 runtime。 |
| 几何 / Volume | ⚠️ `tumour_geometry_dynamic_block` 仍是 logistic + 滤波组合，`tumour_volume_l`、`tcell_density_per_ul` RMSE 大。 | 把几何 ODE 写成 `geometry_whitebox.py`，与 T cell ODE 共用 stiff solver。 |
| Exporter 语义 | ⚠️ 已用 `scripts/dump_snapshot_semantics.py --keyword T1 --keyword V_T` 验证 Reaction 5–27 / volume 规则存在，但遇到缺项仍需补 MATLAB exporter。 | 保持快照完整；几何白盒动工前再覆核一次。 |
| CLI/诊断 | ✅ `scripts/validate_surrogate.py --emit-diagnostics --dump-flat-debug` 可提供所有上下文信号。 | 白盒模块接入后继续用它衡量 A‑series 数值门。 |

---

## 2. 关键资产 / 参考

| 资产 | 说明 |
| --- | --- |
| `matlab/scripts/dev_pd1_training_probe.m` / `scripts/dev_pd1_probe_diff.py` | PD‑1 probe/diff 模板；T cell 版已经落地为 `dev_tcell_probe.m` / `dev_tcell_probe_diff.py`，几何后续照此实现。 |
| `scripts/export_pd1_clean_training.py` | 证明“clean ODE 训练集”可行，后续可仿造出 `export_tcell_clean_training.py` / `export_geometry_clean_training.py`。 |
| `docs/tcell_driver_plan.md` | 描述 surrogate vs MATLAB 的 T cell 行为差异，为白盒化提供背景。 |
| `scripts/dump_snapshot_semantics.py` + `matlab/scripts/show_alignment_drivers.m` | 用来确认 Reaction 14–27、几何规则是否被 exporter 完整写到 snapshot。 |
| `src/offline/stiff_ode.py` | 统一的 solver glue；T cell/geometry 白盒直接复用。 |

---

## 3. 历史坑 / 教训

1. **legacy parquet ≠ ODE**  
   - PD‑1 阶段我们已经确认：`artifacts/training/pd1_whitebox_training.parquet` 中的 state/`pd1_inhibition` 是经过 MATLAB 端的后处理（clip/filter）；直接拿来作为 RMSE 硬目标会逼着白盒去复刻旧 hack。  
   - 对 T cell/几何必须一开始就明确：只用 probe CSV 和 clean exporter 来做物理对齐，legacy CSV 仅用于回归或 sanity check。
2. **单位 / 深度转换**  
   - PD‑1 曾经把 `kon` 多乘了 `86400` / `1e3`，导致 flux 相差 10⁵。T cell/volume 的 `vol_cell`, `k_T1_*`, `geom_*` 同样需要对照 `equations.txt` 和 `parameters.csv` 仔细转换，避免重复踩坑。
3. **Exporter 漏项**  
   - 若 snapshot 缺某条 ODE（例如 `V_T` 只是 repeated assignment），Python 白盒就算写出来也没法初始化/写回。必须在写代码前先跑 `scripts/dump_snapshot_semantics.py` 确认所有反应、规则、事件都能从 snapshot 获得；若没有，先补 MATLAB exporter 再开工。
4. **求解器参数复用**  
   - stiff solver 的 `max_step_days/rtol/atol` 已经在 PD‑1 上验证，不要在 T cell/geometry 模块里另起炉灶或擅自调小步长，否则 A-series CLI 运行会再次出现“多重 solver 配置”导致的难以 debug 行为。
5. **诊断信号**  
   - 过去 PD‑1 调试时忘记在 context 写回中间量，导致 `--dump-flat-debug` 看不到 state；本次 T cell/geometry 白盒要确保 `writeback` 中写出全部 observable（如 `tcell_density_per_ul`, `geom_volume_*`），以便 CLI/metrics 使用。

---

## 3. 行动计划

### Step A — 梳理方程与参数
1. 用 `scripts/dump_snapshot_semantics.py artifacts/matlab_frozen_model/example1 --keyword T1 --keyword nT1 --keyword V_T` 提取 T cell、几何相关的 ODE / 规则 / 事件。
2. 把所有常数（`k_T1_act`, `k_T1_pro`, `q_T1_*`, `geom_growth_per_day`, `vol_cell` 等）收集进新的 dataclass（`TCellParams`, `GeometryParams`），并确认 snapshot 中存在对应条目；若没有，补 exporter。
3. 明确 observable 定义：`tcell_density_per_ul`, `tumour_volume_l`, 以及任何组合信号（伪进展指数、死亡体积）都要写明公式。

### Step B — MATLAB probe
1. ✅ 已完成：`matlab/scripts/dev_tcell_probe.m` 现已导出 `nT1/aT1/T1/T_exh`、`H_PD1_C1`、`H_mAPC`、`H_P1` 以及所有 finite diff；样例输出见 `artifacts/dev/tcell_probe_pd1_train_0004.csv`。
2. ✅ `dev_geometry_probe.m` 可生成 `tumour_volume_l` + live/dead/T cell 体积；配套 `scripts/dev_geometry_probe_diff.py` 在 `pd1_train_0002/0004` 上已验证几何白盒与 MATLAB 对齐。（`pd1_train_0001` 仍受 MATLAB solver 容忍度限制，暂跳过。）

### Step C — Python diff + 白盒模块
1. ✅ 已完成初版：`src/offline/modules/tcell_whitebox.py` 已加载 snapshot 参数并接入 stiff solver，但 RHS 仍需按 probe diff 校正，单测框架待补。
2. ✅ `scripts/dev_tcell_probe_diff.py` 现支持 `--track-derivs` 与 flux 记录；`tcell_whitebox.py` 也已匹配 MATLAB 的 `(Tumor_T1+Tumor_T0)/V_T` 观测，`pd1_train_0004` diff 显示 `RMSE_density≈1.5e-7`、`RMSE_dTumor_T1_dt≈2e-11`。
3. ✅ 几何模块同样具备 probe/diff (`geometry_whitebox.py` + `scripts/dev_geometry_probe_diff.py`)，`pd1_train_0002/0004` 的 `tumour_volume_l` 与 live/dead cell 轨迹皆达机器精度。

### Step D — Runtime 接入
1. ✅ `alignment_driver_block` 已可在 `alignment_mode=1` 下启用 `tcell_alignment_use_whitebox` / `geometry_alignment_use_whitebox`；`alignment_mode=2` 会同时启用 PD‑1 白盒，但当前在 A1 上触发 `NumericsError`（PD‑1 ODE “Required step size is less than spacing between numbers”）。已计划的修复步骤：
   - 在 `alignment_driver_block` 中缓存 `last_pd1_update_time`，若 `dt ≤ min_dt_pd1`（例如 1e‑6 day）则跳过 `pd1_model.step`，避免在同一时刻被多次调用。
   - 构造 `PD1WhiteboxModel` 时传入宽松一点的 SolverConfig（专属 `rtol/atol/max_step`），并将每段 `dt` 限制为 `min(max_internal_step, dt/4)`。
   - 在 `stiff_ode.solve_stiff_ivp` / `integrate_local_system` 中增加可选 debug hook，记录 `(t0, t1, nfev, njev, status, message)`，以便 future diagnose。
   - 若上述仍有 sporadic failure，捕获 `NumericsError` 时落回解析 steady state 或灰盒 Hill 近似，以 warning 方式继续 runtime。
2. ⏳ CLI 验证：`python -m scripts.validate_surrogate --scenarios A1 --ic-mode snapshot --emit-diagnostics --dump-flat-debug 5 --module-block alignment_driver_block --param-override alignment_mode=1 --param-override tcell_alignment_use_whitebox=1 --param-override geometry_alignment_use_whitebox=1` 已跑通，但因 PD‑1 仍是灰盒滤波，`pd1_occupancy` rel_err≈1e2。待上述 PD‑1 调用/容差修复完成后，再切回 `alignment_mode=2` 重跑并观察 `pd1_occupancy`/`tcell_density_per_ul`/`tumour_volume_l` 三个观测。
3. 若 solver 行为稳定，再扩展到 A2–A6，记录 `artifacts/validation/metrics.csv` 中的变化。
4. 🔄 2025‑11‑19 更新：`alignment_driver_block` 现已直接调用 `_project_pd1_to_synapse` 聚合 `V_C/V_T/V_P/V_LN` 的 nivolumab 浓度，新的调试字段 `pd1_alignment_concentration_pk/projection_molar/projection_surface`、`pd1_whitebox_blocked_fraction` 已在 `--dump-flat-debug` 中输出。A1 运行显示投影浓度约 2.2×10⁻⁷ M，与 MATLAB `drug_tumor_molar` 同级；然而白盒 `syn_pd1_pdl1/syn_pd1_pdl2` 仍保持 0，`H_PD1_C1` 也一直为 0。下一步需要开启 `pd1_alignment_debug_solver=1` 逐段记录 `(t0,t1)` 和 RHS，以确认是否有参数单位或 state writeback 漏洞导致复合物始终无法积累。
5. 🔄 2025‑11‑19（二次更新）：为了解决 repeated assignment 把 `time_days` 反复打回 0 的问题，在 runtime 里新增 monotonic “effective time” 和 `pd1_alignment_pending_dt/step_dt/step_count` 诊断字段，并确保 `syn_T1_C1.*` 与 `syn_pd1_*` 同步写回。现在 `pd1_alignment_step_dt` 只要累计到 `max(pd1_alignment_min_dt_days, 1e-6)` 就会触发 BDF 步进，CLI 日志能看到首次 `should=True` 的调用，但仍然在 `dt≈1e-7` 附近触发 `NumericsError: Required step size is less than spacing between numbers`。下一步需要：
    - 把 `pd1_alignment_min_dt_days` 或内部 pending gate 提到 1e-4 左右，避免 BDF 在 machine‑eps 级别挣扎；
    - 或者在 `pd1_whitebox_model.step` 内部添加兜底，若 span 太小就直接重用上一状态，防止 driver 卡死。
    - 一旦 solver 稳定，再重跑 `--dump-flat-debug` 验证 `syn_pd1_pdl1` 是否离开 0。
6. ✅ 2025‑11‑20 更新：在不改 `switches.py` 主逻辑的前提下，通过 CLI 覆盖 `pd1_alignment_max_step_days=0.01` 与 `pd1_alignment_solver_atol=1e-10`，并启用 `alignment_mode=2` 让三大白盒一起运行，A1 终于可以稳定集齐 `pd1_alignment_step_status=1`（参见 `artifacts/validation/A1_flat_debug_20251119T154410.csv`，共 49 次成功步进、0 次 NumericsError）。运行命令：

```bash
python -m scripts.validate_surrogate --scenarios A1 --ic-mode snapshot \
  --module-block alignment_driver_block --emit-diagnostics --dump-flat-debug 50 \
  --param-override alignment_mode=2 \
  --param-override pd1_alignment_use_whitebox=1 \
  --param-override tcell_alignment_use_whitebox=1 \
  --param-override geometry_alignment_use_whitebox=1 \
  --param-override pd1_alignment_max_step_days=0.01 \
  --param-override pd1_alignment_solver_atol=1e-10
```

上述 run 说明 pending gate=0.5 天时每一步都会被拆分成 solver 内部 0.01 天的小步，`syn_pd1_pdl1≈2.7e-6` 很快“破零”，也为后续 A2–A6 的 numerics 验证提供了可复现模板。

### Step E — 数据/文档
1. 决定是否需要干净的 T cell/几何训练集：如果要拟合，就仿照 PD‑1 写 clean exporter 并生成 parquet。
2. 更新 `docs/new_alignment_plan.md` / `docs/project_handoff.md`：标记 PD‑1 完成、T cell/几何正在白盒化；把 “training parquet RMSE <1e‑2” 归档到 optional/backlog。
3. 维护 `docs/tcell_driver_plan.md`：记录 probe diff 结果、参数调优记录。

---

## 4. MATLAB / 生物团队需确认

1. **参数来源**：`k_T1_act`, `k_T1_pro`, `geom_*` 是否已有权威值？请在 snapshot `parameters.csv` 中补齐，避免 Python 端 hardcode。
2. **Observable 计算**：`tcell_density_per_ul`, `tumour_volume_l` 在 MATLAB 端是否经过额外几何缩放或滤波？请提供公式以免白盒输出错位。
3. **Exporter 补丁**：若 `V_T` 或 `T` 相关的 ODE/规则未导出，请提供脚本改动，确保 snapshot → Python 的方程闭环。

---

## 5. 参考命令

```bash
# Probe + diff（PD‑1 示例，T cell/几何照抄）
/Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab -batch \
  "cd('/Volumes/AlbertSSD/Program/new/qspiopy'); addpath(fullfile(pwd,'matlab','scripts')); \
   dev_pd1_training_probe('pd1_train_0004');"
python scripts/dev_pd1_probe_diff.py artifacts/dev/pd1_training_probe_pd1_train_0004.csv

# T cell probe + diff（已落地）
/Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab -batch \
  "cd('/Volumes/AlbertSSD/Program/new/qspiopy'); addpath(fullfile(pwd,'matlab','scripts')); \
   dev_tcell_probe('pd1_train_0004');"
python scripts/dev_tcell_probe_diff.py artifacts/dev/tcell_probe_pd1_train_0004.csv

# Alignment driver 快速验证
python -m scripts.validate_surrogate --scenarios A1 --ic-mode snapshot \
  --module-block alignment_driver_block \
  --emit-diagnostics --dump-flat-debug 5 --max-rel-err 1e12
```

遵循以上步骤，T cell 与动态几何模块都会像 PD‑1 一样进入 stiff solver 白盒体系；届时 `alignment_driver_block` 只需做模块拼装，灰箱调参历史问题也将随之消失。
