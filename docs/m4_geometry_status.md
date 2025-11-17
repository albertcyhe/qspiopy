# M4 动态几何对齐 — 快速入门（2025‑11‑10）

本文是当前 M4（动态几何 + 医疗 PK/PD 调度）状态的自包含说明，适合刚接手的同学快速了解「已经完成的」「如何复现实验」「尚待解决的问题」以及「下一步优先级」。

---

## 1. 当前进度概览

| 模块 | 状态 | 备注 |
| --- | --- | --- |
| **几何实时模块** | ✅ | `tumour_geometry_dynamic_block` 已更新为直接使用 `context['tumour_volume_l']` / `context['V_T']` 作为基线，再叠加 live/T-cell/dead 滤波体积；自动识别 `vol_cell` 单位，并在 `--dump-flat-debug` 下输出 `geom_*` 诊断键。 |
| **模块注入策略** | ✅ | 自 2025‑11‑13 起，A 系列 CLIs 默认启用 `alignment_driver_block`（白盒 PK / PD‑1 / 几何 ODE）。旧的 `pd1_bridge_block`、`tumour_geometry_dynamic_block` 仍可通过 `--module-block` 明确指定，用于 debug 与回归。 |
| **CLI 诊断** | ✅ | `scripts/validate_surrogate.py` / `--emit-diagnostics` / `--dump-flat-debug` 已覆盖 snapshot & target-volume 两种初始化模式。 |
| **对白盒 PD‑1 的依赖** | 🔄 | Alignment driver 的 PD‑1 子模块已改为显式 kon/koff/k_int ODE，但默认参数与训练曲线差异仍大；需要结合 PD‑1 白盒的最新进展一起调参。 |

---

## 2. 关键资产

| Path / 命令 | 用途 |
| --- | --- |
| `src/offline/modules/switches.py` (`tumour_geometry_dynamic_block`) | 减缓肿瘤体积、伪进展滤波、T cell density 计算。 |
| `scripts/validate_surrogate.py` | 主 CLI。`--module-block ...` 控制是否启用旧几何模块，`--emit-diagnostics` + `--dump-flat-debug 5` 输出对齐所需数据。 |
| `artifacts/matlab_frozen_model/example1/parameters.csv` | 参考的 `k_cell_clear`, `vol_cell`, `vol_Tcell` 等几何常数。 |
| `docs/alignment_tuning_plan.md` | 更上层的调参路线（与 PD‑1 / 几何位置关系）。 |

---

## 3. 如何复现当前结果

### 3.1 Snapshot 路径

```
python -m scripts.validate_surrogate --scenarios A1 --ic-mode snapshot \
  --emit-diagnostics --numeric-gates \
  --module-block pd1_bridge_block \
  --module-block pd1_occupancy_filter_block \
  --module-block tumour_geometry_dynamic_block \
  --dump-flat-debug 5
```

关键指标（目前仍不合格）：

| Observable | rel‑L2 | maxRE | 备注 |
| --- | --- | --- | --- |
| `tumour_volume_l` | ≈ 7.9e‑1 | ≈ 1.27 | 基线 1.4e‑2 L，但下行趋势与 MATLAB 差异大。 |
| `tcell_density_per_ul` | ≈ 6.1 | ≈ 3.2e+1 | 需要使用 `V_T.T*` 聚合或额外驱动。 |
| `pd1_occupancy` | ≈ 1.0 | ≈ 5.6e3 | `pd1_occupancy_filter_block` 现为「当前 Hill 输出」，尚未调成 MATLAB 的延迟滤波。 |

### 3.2 Target-volume 初始化

```
python -m scripts.validate_surrogate --scenarios A1 --ic-mode target_volume \
  --ic-target-diam-cm 0.05 --ic-max-days 400 --ic-max-wall-seconds 60 \
  --emit-diagnostics --numeric-gates \
  --module-block pd1_bridge_block \
  --module-block tumour_geometry_dynamic_block \
  --dump-flat-debug 5
```

- 0.5 cm / 150 d 的默认 IC 无法在时限内收敛；改为 0.05 cm / 400 d 可以完成，但仿真输出仍近似水平（`tumour_volume_l ≈ 1.78e-4 L`，`tcell_density_per_ul ≈ 4e5`），数值门依旧超标。

---

## 4. 当前阻塞

1. **PD‑1 占有率**  
   alignment driver 的 PD‑1 ODE 未能复现 MATLAB 的缓慢爬升（rel‑L2≈O(1)）。需要配合 PD‑1 白盒 fitter 调整 `kon/koff` 缩放、`PD1_50`、滤波时间常数。
2. **肿瘤体积 / T cell density**  
   几何 follower 的 logistic 参数沿用默认值，导致 volume / density 在数天内剧烈震荡。需要重新估计 `geom_growth_per_day`, `geom_kill_per_cell_per_day`, `geom_volume_cap_l` 等参数，或直接从 MATLAB 导出等效 ODE。
3. **Exporter 语义不完整**  
   若 MATLAB snapshot 中 `V_T` 通过 repeated assignments 直接写入（而非 ODE），则 Python 的白盒路径无法独立演化，需要补全导出的方程或额外 metadata。

---

## 5. 需要 MATLAB 同步的信息

1. **几何参数的真实值**  
   `k_cell_clear`, `vol_cell`, `vol_Tcell`, `geom_*` 等是否有明确来源？请在 MATLAB 导出前写入 snapshot（`parameters.csv`），让 Python 自动拾取。
2. **PD‑1 占有率是否有额外延迟/事件**  
   如果 `H_PD1` 的滤波依赖其它变量（例如 tumour volume, surface area scaling），请提供公式/事件描述，便于 Python 模块复现。
3. **初始条件**  
   如果 MATLAB 参考曲线的体积/密度不是从 `V_T` 直接读取，需要明确“真实输出”是哪个变量，以便 Python 端取用一致的 observable。

---

## 6. 下一步（按优先级）

1. **调 `pd1_occupancy_filter_block`**  
   - 读取 MATLAB 的 `H_PD1` 曲线，拟合滤波参数（`tau`, `delay`, `PD1_50`）。  
   - 目标：`pd1_train_0004`/`0582` 的 `H_RMSE < 1e-2`（与 PD‑1 白盒同一验证逻辑）。
2. **重新估计几何参数**  
   - 从 MATLAB 模型或文献获得 `geom_*` 建议值，并更新 `parameters/example1_parameters.json`。  
   - 通过 `scripts/validate_surrogate.py --emit-diagnostics` 验证 volume / density RMSE 是否下降。
3. **Exporter 补全**  
   - 确认 MATLAB ODE 是否完整导出；若没有，则增补 exporter 使 Python 可独立运行。  
   - 根据结果调整 `alignment_driver_block` 逻辑，并决定何时完全切换到纯 snapshot 模式。

完成上述三项后，再回到 `scripts/dev_pd1_driver_compare.py` / `scripts.validate_surrogate --scenarios A1…A6`，观察全部 observable 是否满足数值门，然后更新 `docs/test_status.md`。
