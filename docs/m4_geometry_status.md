# M4 动态几何校准现状（2025-11-10）

## 已完成的 Python 侧改动

1. **动态几何模块**：`tumour_geometry_dynamic_block` 负责聚合多克隆体积并对死亡细胞做一阶低通滤波（伪进展）。最新实现会：
   - 自动识别 `vol_cell` 是否已是升单位（≤1e-6）或仍为 µm³，并据此换算。
   - 在运行期使用 `context['tumour_volume_l']` / `context['V_T']` 作为基线体积，再叠加 live/T-cell/dead 滤波体积，避免重新计算后得到 1e-17 L 级别的数值。
   - 暴露 `geom_*` 诊断键供 `--dump-flat-debug` 查看。
2. **模块导出与 CLI**：`pd1_bridge_block` + `tumour_geometry_dynamic_block` 原本是 A1 默认模块；自 2025‑11‑13 起，A 系列默认改为只启用 `alignment_driver_block`（更白盒的 PK/PD‑1/体积 ODE）。旧模块仍留作 debug/回归引用。

## Snapshot 路径（ic_mode=snapshot）结果

命令：
```
python -m scripts.validate_surrogate --scenarios A1 --ic-mode snapshot \
  --emit-diagnostics --numeric-gates --module-block pd1_bridge_block \
  --module-block pd1_occupancy_filter_block \
  --module-block tumour_geometry_dynamic_block --dump-flat-debug 5
```

关键观测：
- `tumour_volume_l`: rel_L2 ≈ **7.9e-1**, maxRE ≈ **1.27** — 基线体积为 1.4e-2 L，但仍缺少 MATLAB 式的缓慢下降。
- `tcell_density_per_ul`: rel_L2 ≈ **6.1**, maxRE ≈ **3.2e+1** — 分子 (T_total) 仍近乎常数，需在灰箱里把 `T_total` 替换为写入的 `V_T.T*` 聚合或另做驱动。
- `pd1_occupancy`: rel_L2 ≈ **1.0**, maxRE ≈ **5.6e3** — 新增的 `pd1_occupancy_filter_block` 已接管输出，但尚未调参，当前等效为“近期 Hill 值”，与参考缓慢爬升差异巨大。

## Target-volume 初始化尝试

命令：
```
python -m scripts.validate_surrogate --scenarios A1 --ic-mode target_volume \
  --ic-target-diam-cm 0.05 --ic-max-days 400 --ic-max-wall-seconds 60 \
  --emit-diagnostics --numeric-gates --module-block pd1_bridge_block \
  --module-block tumour_geometry_dynamic_block --dump-flat-debug 5
```

现象：
- `generate_initial_conditions` 无法在 150 d 内达到默认 0.5 cm（最后只有 0.064 cm），改为 0.05 cm / 400 d 能完成，但仿真输出仍然近似水平（`tumour_volume_l ≈ 1.78e-4 L`，`tcell_density_per_ul ≈ 4e5`），数值门依旧超标。

## 当前瓶颈（更新）

1. **PD‑1 占有率**：`alignment_driver_block` 已重写为显式 kon/koff/k_int ODE，但默认参数仍无法复现 MATLAB 的缓慢爬升（rel‑L2≈O(1)）。需要借助 `scripts/fit_observables.py` 拟合新参数。
2. **肿瘤体积 / tcell 密度**：alignment driver 的 logistic follower 按默认值会在数天内大幅波动。必须调 `geom_growth_per_day`、`geom_kill_per_cell_per_day`、`geom_volume_cap_l` 等新引入的白盒参数。
3. **Exporter 语义**：若 MATLAB 端未导出完整的几何/PD‑1 ODE（例如直接在 repeated assignment 里写 `V_T`），需要修改 exporter 并重新冻结 snapshot，确保 Python 白盒路径能独立运行。

## 需要 MATLAB 专家协助的问题

1. **几何参数的权威数值**：SimBiology 模型中是否有与 `k_cell_clear`, `vol_cell`, `vol_Tcell` 搭配使用的常数，可直接转化为上述 `geometry_*` 参数？若有，请在 MATLAB 导出前写入、以便 Python 自动拾取。
2. **PD‑1 事件**：MATLAB 端的 PD‑1 占有率是否依赖其它触发（例如 tumour volume, cell surface area scaling）？若是，需要确认导出的 repeated assignment/模块在 Python 端是否完整执行。
3. **初始条件**：MATLAB 参考文件中的 `tumour_volume_l` 波形看似非零起点且有缓慢变化；请确认在 MATLAB 中取用的是 `V_T`（或其它变量），并说明是否存在事件/赋值在 t=0 即重新设定体积。

## 建议的后续步骤（待 MATLAB 回答）

1. 获取/确认几何参数在 MATLAB 端的真实值，并写入 snapshot（`artifacts/matlab_frozen_model/<snapshot>/parameters.csv`）。
2. 若 MATLAB 端对 PD-1 占有率还有额外的延迟/平滑处理，请描述相关公式，以便在 Python 模块中复现。
3. 如果 snapshot ODE 无法在 150–400 d 内长成 0.5 cm，需要 MATLAB 给出推荐的 target diameter / integration horizon，或在 snapshot 中直接提供“已长成”状态。

> 注：一旦 exporter/参数同步完成，`alignment_driver_block` 将退化为可选 debug 工具，A 系列可以回到“纯 snapshot”语义。

附：相关文件与命令
- 动态几何模块：`src/offline/modules/switches.py`（`tumour_geometry_dynamic_block`）
- CLI 入口：`scripts/validate_surrogate.py`
- 参数检查：`parameters/example1_parameters_inspect.csv`
- 诊断命令：见上述两条 CLI 示例。
