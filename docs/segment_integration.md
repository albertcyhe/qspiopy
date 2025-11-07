# QSP-IO 分段积分与事件处理方案

本文档说明了目前在 Python 端对 SimBiology 快照进行积分时采用的“分段积分 + 事件重演”策略，以及后续扩展方向。

## 设计目标

1. **显式处理离散事件**：包括 bolus 给药、快照自带的 SimBiology Event（trigger/delay）以及延迟触发的作业。
2. **避免在断点处硬踩求解器**：在事件/给药时刻前结束连续积分，事件执行后再从 `t_event + ε` 起步，防止出现 “required step size …” 的拒步。
3. **保持记录语义**：对事件发生前后的状态都应存盘，以便与 MATLAB/SimBiology 输出一一对齐。

## 核心做法

1. **构造断点序列**  
   将 0、所有给药时间、仿真结束时间 `days`（去重、排序）组成 `segment_points`。连续积分只在 `[segment_i, segment_{i+1})` 区间内进行。

2. **分段调用 `solve_ivp`**  
   - 每段的上界初始为下一断点，如存在已排队的延迟事件，则取二者较小者。  
   - SciPy 求解器始终工作在连续区间内，不携带 `t_eval`，依靠 dense output + 事后采样。
   - 若 `segment_end` 与断点重合，则在数值上后推一 ULP（`np.nextafter(target, -inf)`），避免求解器需在断点评估。

3. **事件检测与重演**  
   - SimBiology Event 通过 `solve_ivp(..., events=...)` 触发；检测到事件后，记录事件前状态 → 执行赋值/排入 delay → 记录事件后状态。  
   - delay 事件（`event.delay_type == "time"`）插入最小堆 `pending_events`，下一轮积分前会把上界截断到最早的 delay。

4. **重复赋值、代数规则与采样**  
   - 每次积分返回的状态都会通过 `reconcile()` 重新评估 repeated assignment 和 algebraic rule，保证派生变量（PD-1 占有率等）与 MATLAB 语义一致。  
   - 采样时刻使用 `_merge_close_times` 去除 1 ULP 内的重复，并在字典键上统一 `_tkey()` 量化，避免浮点误差造成重复或遗漏。

5. **剂量执行**  
   - 到达断点后，对所有排在该时刻的 `ScheduledDose` 执行 `apply_dose()`，随即记录事件后状态，并用 `np.nextafter()` 将时间推进至断点右侧作为下一段起点。

## 仍需改进

1. **与求解器更紧密的“微步”策略**：当前只在断点处回退一 ULP，如仍由模型事件触发拒步，需要再加入“事件触发 → 直接跳入断点 → 事件后起步”的容错路径。
2. **通用的分段积分模块化**：后续可以把现有逻辑抽成 `segment_integrator.py`（可供离线 runner、CI、以及其他 CLI 复用），并将求解器配置（method/rtol/atol/jacobian）包装成不可变结构做 provenance。
3. **单元 / 组件测试**：需要为“多事件同刻”、“delay 事件”、“连续 dose”构造最小快照并在 CI 中回归，确保新的分段策略不会回退。
4. **进一步的性能优化**：在大规模场景（人群、灵敏度分析）中，可以通过传递稀疏 Jacobian 或切换到 `numba`/`jax` 生成的 RHS 来降低分段开销。

## 与 SimBiology 的语义对齐

| SimBiology 行为               | Python 实现                                                         |
|------------------------------|---------------------------------------------------------------------|
| 事件触发（上升沿）          | `solve_ivp(events=...)` 检测，内部将 trigger 转换成零残函数         |
| 事件执行时刻微偏移          | 断点积分到 `target^-`，再手动在 `target` 上执行赋值并推进 `+ε`      |
| 事件前/后状态各一行         | `_record_solution_samples` + `samples[_tkey(event_time)] = state`    |
| delay 事件                  | `pending_events` 最小堆，按触发时间二分积分上界                    |
| repeated assignments / rules | `reconcile()` 在每次积分 / 事件后重评，并同步回 state/context       |

## 文件位置

- 分段积分逻辑：`src/offline/simulation.py`（主循环中按断点推进的部分）
- 快照元数据与 repeated assignment：`src/offline/snapshot.py`
- 后续若要提取为独立模块，可参考 `docs/software_architecture.md` 中的“策略点”约定。

---

如需复现此设计，可在 `simulation.py` 中查看 `segment_points`/`pending_events` 以及与 `solve_ivp` 的交互；若要在 CLI / CI 中验证，可运行：

```bash
python -m scripts.validate_surrogate --scenarios A1 --dump-t0 --numeric-gates
```

它会在分段积分基础上执行 t=0 快检、剂量审计以及数值门禁。
