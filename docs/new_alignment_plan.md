之前**段式积分器**已经把“卡在剂量断点”的稳定性问题解决了，现在，系统性地把 **SimBiology ↔ Python** 的语义/单位/模块差异一口气对齐起来。下面给出一份**工程化的改造方案**（包含分阶段目标、文件级改动点、关键函数签名、示例代码片段、测试与验收要点）。你可以按阶段落地，每个阶段都能独立带来实打实的对齐收益。

---

## 总览（里程碑与进度）

| 里程碑 | 目标 | 当前状态 |
| --- | --- | --- |
| **M1 事件语义** | 状态触发 + 同刻顺序 + pre/post 双记录 + ε‑bump + first_step + 事件后 reconcile | ✅ 已实现：schedule 清洁、聚合同刻剂量、延迟事件队列、trigger specs、metadata (`phase_code` 等) |
| **M2 单位与参数** | `units.py` 统一换算（时间/体积/浓度/速率/剂量），参数派生 | ✅ 已完成：所有时间/体积/流量/kon/koff 路径统一到 day/L/M；`normalise_dose_to_species()` 驱动 `apply_dose()`；ParameterGraph 派生值写入 `unit_normalisation_map` 并可由 `scripts/print_units_table.py` 审计；新增 `tests/test_units.py`、`tests/test_param_graph.py`。**遗留风险**：2D kon 仍依赖 legacy 常数（待几何参数化）；快照若缺药物 MW 则会在 mg 剂量路径上硬 fail；A1 数值门虽然跑通流程但 tumour/occupancy/tcell_density 仍 ❌（语义问题挪至 M3/M4 解决）。 |
| **M3 初始化与模块化** | 目标体积初始条件、模块化加载 | 🟡 进行中：`simulate_frozen_model` 默认 `ic_mode="target_volume"`（0.5 cm、150 d、`reset_policy="cancer_only"`；不支持的快照自动退回 snapshot 初值），CLI 暴露 `--ic-mode/--ic-target-diam-cm/--ic-max-days/--ic-max-wall-seconds`；state-trigger 事件已加入去抖/迟滞状态机（armed/disarmed + 不应期 + 回到 false 再触发），并配套 `tests/test_state_trigger_hysteresis.py` 与 `tests/test_events.py -k event_suite` 覆盖。**残留**：尚未在完整 snapshot 场景（如 example2/A1）上复核；确认稳定后再把 `ic_mode=snapshot` 设为默认并重跑数值门。 |
| **M4 多克隆与动态体积** | 体积/伪进展输出 & 克隆竞争 | ⏳ 未开始 |
| **M5 验收/CI** | 组件测试 + 数值门绿灯 + CI | ⏳ 进行中（validate_surrogate 现已稳定，但 A1 数值门仍未过） |

- [x] `validate_surrogate` 默认关闭性能基准（`--benchmark-replicates=0`）。
- [ ] 数值门：A1 仍超标（tumour_volume/pd1_occupancy/tcell_density；PK 尾部），待完成 M3/M4 语义梳理后重跑 `--numeric-gates`。

**最新实测（2025-11-07）**

- MATLAB 参考已用 `/Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab` 重新生成（`python -m scripts.run_alignment_suite --scenarios A1 --output artifacts/validation`）。
- `python -m scripts.validate_surrogate --scenarios A1 --dump-t0 --numeric-gates`：当前仍以 `ic_mode=snapshot` 运行以便直接对照 frozen snapshot（默认 CLI 已改为 `ic_mode=target_volume`；snapshot 路径需等迟滞补完后再启用），tumour_volume_l rel_L2≈1.4e-1、pd1_occupancy rel_L2=1.0、tcell_density rel_L2≈0.79。待状态触发事件的迟滞语义在轻量诊断场景验证通过后，再用 `--ic-mode snapshot` 复测；在此之前，验证/runner 默认使用 `--ic-mode target_volume --ic-target-diam-cm 0.5 --ic-max-days 150 --ic-reset-policy cancer_only`。

**M3 剩余重点**

1. 状态触发事件迟滞/防抖：触发一次后必须“退回安全区”再触发（包含 ε‑bump、armed state、时间抑制窗）。
2. 轻量诊断场景：构造最小 ODE + 状态事件测试，确保迟滞逻辑不会在 t≈0 重复触发。
3. Snapshot IC 复用：上述两项完成后，再把 `ic_mode=snapshot` 作为默认，并重新跑 A/B 系列。

---

## M1：事件语义 100% 对齐

### 1) `src/offline/segment_integrator.py`（**已有**）增强

**目标**：不只是“定时剂量”，还要**状态触发型事件**、并发事件顺序、同刻 **pre/post** 双记录、ε‑bump 与受控 `first_step`、事件后 reconcile。

**新增/调整 API**：

```python
@dataclass(frozen=True)
class ScheduledDiscontinuity:
    # 统一表示“剂量 / 立即事件 / 延迟事件触发”
    time: float
    kind: Literal["dose", "event_immediate", "event_delayed"]
    priority: int           # 使用 event_index 或 dose_index 保序
    payload: Any            # DoseEntry 或 EventEntry 及其赋值列表等

@dataclass(frozen=True)
class TriggerEventSpec:
    # 状态触发事件，用 solve_ivp events 机制找零点（上升沿）
    entry: EventEntry
    direction: float        # from snapshot.EventEntry.direction
    fn: Callable[[float, np.ndarray], float]  # returns signed distance

def run_segmented_integration(
    rhs: Callable[[float, np.ndarray], np.ndarray],
    y0: np.ndarray,
    t_span: tuple[float, float],
    sample_times: np.ndarray,
    *,
    # ① 固定时间断点（剂量/已调度事件/延迟事件）
    schedule: Sequence[ScheduledDiscontinuity],
    # ② 状态触发事件
    triggers: Sequence[TriggerEventSpec] = (),
    # ③ Jacobian / 稀疏模式
    method: str = "BDF",
    rtol: float = 1e-7,
    atol: float = 1e-10,
    max_step: float = np.inf,
    jac_sparsity: Optional[np.ndarray] = None,
    # ④ 语义对齐钩子
    apply_discontinuity: Callable[[ScheduledDiscontinuity, np.ndarray], np.ndarray],
    reconcile: Callable[[np.ndarray], None],   # repeated/algebraic rules
    record: Callable[[float, np.ndarray, Literal["pre","post","cont"], Optional[int], Optional[str]], None],
    # ⑤ 其它策略
    eps: float = np.finfo(float).eps,
) -> None: ...
```

**核心实现要点**：

* 事件**触发**：将 `EventEntry.trigger_compiled` 转换为 `TriggerEventSpec.fn`；设置 `fn.direction = entry.direction`，`fn.terminal = False`。把 `triggers` 传给 `solve_ivp(events=...)`。
* **同刻并发事件顺序**：按 `(time, priority)` 排序，`priority = event_index`（或 `dose_index`）。同刻事件依次执行，每次执行后立刻 `reconcile`，可连锁触发下一事件（SimBiology 语义）。
* **pre/post 双记录**：每次断点执行：

  * 事件/剂量**前**：`record(t, y_pre, "pre", event_index, entry.name)`
  * 执行离散赋值 → `reconcile` → `y_post`
  * **后**：`record(t, y_post, "post", ...)`
* **ε‑bump**：每次断点后都从 `t_next = np.nextafter(t, +∞)` 作为下一段的起点；并将 `first_step` 设为 `min(0.5*max_step, max(10*eps, 0.01*(t_target-t_next)))`。
* **挂死保护**：若 SciPy 报 *“required step size is less than spacing between numbers”* 或 *“ts must be strictly increasing”*，直接：

  1. 记录边界 pre/post；
  2. bump 到 `t + ε`；
  3. 重启一段。如果连续三次遇到相同的边界错误，抛 `NumericsError`，并将当前 `schedule` 条目写入错误消息（避免静默循环）。

> 你已经做了 ε‑bump/受控 first_step/样本密集写入，我建议把 **pre/post 双记录** 改为**由 integrator 统一完成**（通过 `record` 回调向 `ScenarioResult.extras` 写 `phase="pre"/"post"`），减少 `simulation.py` 里重复代码。

---

### 2) `src/offline/snapshot.py`（**补齐**）

**新增**：

```python
class FrozenModel:
    ...
    def build_trigger_specs(self) -> list[TriggerEventSpec]:
        specs = []
        for entry in self.events:
            # 只为“状态触发”的事件建立 trigger（时间触发不需要）
            if entry.trigger_expression and "time" not in entry.trigger_expression.lower():
                def make_fn(e: EventEntry):
                    def fn(t, y):
                        # 仿 Echt：从 y -> context -> repeated/algebraic -> 评估 e.trigger_compiled
                        ctx = self.build_context_from_state(y.copy())
                        self.evaluate_repeated_assignments(ctx)
                        self.apply_algebraic_rules(ctx, y, mutate=False)
                        return e.trigger_compiled.evaluate(ctx)
                    fn.direction = e.direction
                    fn.terminal = False
                    return fn
                specs.append(TriggerEventSpec(entry=entry, direction=entry.direction, fn=make_fn(entry)))
        return specs
```

* 仍保留你已加的 `jacobian_sparsity()`。
* 若 `EventEntry.delay_type == "time"` 且 `delay > 0`，触发后在 integrator 里放入 **延迟事件**（`ScheduledDiscontinuity(kind="event_delayed")`），同刻**按 index 顺序**执行。

---

### 3) `src/offline/entities.py`

**目标**：在结果里**显式区分事件前后**，并与“连续采样”统一。

**做法**：不要改主字段（兼容现有消费者）。**通过 `extras` 增列**：

* `phase`: 0=cont（连续）、1=pre、2=post
* `discontinuity_type`: 0=none、1=dose、2=event_immediate、3=event_delayed
* `event_index`, `event_name`, `dose_name`, `target`（可为空）
* `time_key`: 用于调试重复时间写入（量化后时间键）

> 你的 `ScenarioResult.to_frame(order="contract")` 会把 `extras` 一并导出，现有消费者无需变更即可享受更丰富的结果语义。

---

### 4) `src/offline/simulation.py`

**职责收敛**：

* 仅负责：加载 snapshot → 准备 schedule（去重/聚合）与 trigger_specs → 生成 `apply_discontinuity/reconcile/record` 回调 → 调用 `run_segmented_integration`。
* **不要**再手写“事件/剂量 pre/post 记录”，交由 integrator 的 `record()` 完成。

**关键改动点**：

* **Schedule 清洁化**（你已做）：
  ① 去重相同 `(time, target, name)` 的剂量；② 去掉重复 `t=0`; ③ fallback 护栏（若已有等效剂量则不补）；④ 把同刻多剂量聚合为**单次**状态更新（integrator 里聚合已经做了）。
* **触发事件**：
  `triggers = model.build_trigger_specs()`；将它传入 integrator。
* **record 回调**：
  把 `phase/discontinuity_type/event_index/dose_name` 等写入 `ScenarioResult.extras`（维护一个 `samples` 字典用时间键 `_tkey(t)` 管理状态快照）。
* **reconcile 回调**：
  就是你已有的逻辑：`build_context_from_state → evaluate_repeated_assignments → apply_algebraic_rules → sync_state_from_context`。
* **jac_pattern**：
  直接 `jacobian_sparsity()` 传给 integrator；Radau/BDF 会显著少报“奇异”警告。

---

## M2：单位与参数（**关键**）

引入 `src/offline/units.py` 作为**单一可信**的换算层，并**重构** `snapshot._convert_parameter_value/_convert_compartment_value` 与 `FrozenModel.apply_dose()` 的单位逻辑**只调用**该模块。

### 1) `units.py`（建议接口）

```python
from __future__ import annotations
import math
from typing import Literal

TimeUnit = Literal["day"]         # 统一到 day
VolUnit  = Literal["L"]           # 统一到 L
ConcUnit = Literal["mol/L"]       # 统一到 M
AmtUnit  = Literal["mol"]

def time_to_day(value: float, unit: str) -> float: ...
def vol_to_litre(value: float, unit: str) -> float: ...
def amount_to_mol(value: float, unit: str, mw_g_per_mol: float | None = None) -> float: ...
def conc_to_M(value: float, unit: str) -> float: ...

# 反应/速率常用换算
def rate_to_day(value: float, unit: str) -> float: ...           # e.g. 1/s → /day
def kon_to_L_per_mol_day(value: float, unit: str) -> float: ...  # e.g. 1/(µM·nm·s) → L/mol/day
def area_to_m2(value: float, unit: str) -> float: ...
def length_to_m(value: float, unit: str) -> float: ...
```

**处理要点**：

* 全部**统一到 day/L/mol/M**。
* 解决 QSP‑IO 中怪单位：

  * `1/(micromolarity*nanometer*second)` → `1/(µM·nm·s)`：需按 `1 µM = 1e-6 mol/L`、`1 nm = 1e-9 m`、`1 s = 1/86400 day` 推导；**你此前在 `_convert_parameter_value` 用了常数 `9.8412890625`，把推导移到 `units.py` 并覆盖更多变体写法**（大小写、符号、单位拼写变体）。
  * 所有 `1/second`, `1/minute`, `1/(M*s)` 等全部走 `rate_to_day/kon_to_L_per_mol_day`。
* **剂量统一**：`apply_dose()` 改用：

  1. `amount_mol = amount_to_mol(dose.amount, dose.amount_units, mw)`（若 `amount_units` 是 mg 需要药物 MW；把抗体 MW 放 `parameters` 或快照 metadata）；
  2. 若目标物种维度是浓度（从 `SpeciesEntry.interpreted_dimension` 判断），按 `delta_conc = amount_mol / compartment_volume_L`；
  3. 否则按 amount（mol）直接加到数量态。

> 将 `snapshot._convert_parameter_value/_convert_compartment_value` 里的零散转换全部改为 `units.py` 调用。把所有**时间相关速率**（CL、Q、γ、k_on/k_off/k_deg…）统一到 `/day`。

### 2) 参数派生与表达式

**目的**：匹配原文 parameter JSON 的 `derived_from`/`expression` 语义。

* 在 `snapshot.py` 读取参数后，增加 `ParameterManager`：

  * 维护 `{name: value, unit_str: str}`；
  * 解析 `rules.csv` 里**非 ODE 的“参数赋值类表达式”**作为派生（比如 `k_clear = CL/Vc`）；
  * 以 DAG 拓扑序求值；
  * 把派生结果写回 `FrozenModel.parameters`；
  * provenance 里新增 `parameter_units_map`。

---

## M3：初始化与模块化

### 1) 初始条件生成（目标直径/体积）

新增 `src/offline/init_conditions.py`：

```python
def solve_init_to_target_volume(
    model: FrozenModel,
    target_diameter_cm: float,
    *,
    max_days: float = 365.0,
    tol: float = 1e-3,
    method: str = "BDF",
) -> np.ndarray:
    """
    从 model.initial_state 出发，在无治疗（therapy='none'）下演化，
    找到使肿瘤直径接近 target 的状态向量，作为新的初始条件返回。
    """
```

* 内部用 `run_segmented_integration` 长时间段（无离散事件），以 **root‑finding** on time 或**二分**找到最接近目标直径的时间点状态。
* 这样可避免从“非物理解”的 snapshot 初始量出发带来的长暂态偏差（匹配 QSP‑IO 初始化策略）。

### 2) 模块化骨架

* 把 cancer / antigen_APC_MHC / Tcell 拆成策略接口，`FrozenModel` 上暴露 `module_flags`（从快照/配置读取），`simulation.py` 依据 `flags` 选择是否加载对应 ODE/事件/剂量清单。
* **短期**可不大动 ODE，只做“可插拔开关”，减少“缺模块”时产生的对齐偏差。

---

## M4：多克隆 & 动态体积/伪进展

**逐步推进**（避免一次性修改过大）：

1. **体积动态化**：
   在 `simulation.py` 输出处不要再用 `ctx.get("V_T")` 的静态量。改为**根据 C/T/dead 动态计算**：
   `V = (α_C * C + α_T * T + α_D * D)*细胞体积`，参数从快照/配置读取；`diameter` 也由该 `V` 反推。把这一逻辑做成 `ExtraOutputs` 插件（你已有 `ContextKeyOutputs` 框架，可以追加 `TumourGeometryOutputs`）。
2. **伪进展**：
   让 dead cell 清除率影响 `V`（先不改 ODE，仅在输出层表达“体积不立即下降”）。
3. **克隆**：
   将 C/T 分别允许多克隆（`C1..Cn`/`T1..Tm`）；在 `evaluate_ode_rhs` 聚合“总载量”约束，或在 repeated assignment 中更新共享资源项。**第一步**先在输出聚合（`C_total/T_total`）对齐指标；第二步再把克隆间竞争参数引入 ODE。

---

## M5：测试与验收（含 A1 门）

### 1) 组件测试（pytest）

* **事件/延迟**：构造 3 个事件：同刻 `event_index` 递增，包含一个 `delay>0`，断言 pre/post 序列与顺序正确；断言 `phase`/`event_index` extras 正确。
* **剂量微段**：同刻多剂量聚合，断言只产生一次不连续（pre/post 2 行），并且 delta 等于量之和。
* **状态触发**：构造 `y` 穿越阈值触发“上升沿”，断言只触发一次（不反复在边界抖动），bump 后能继续。
* **单位**：对 `units.py` 的每个变体单位写参数化测试（s/min/h/day, µL/mL/L, µM/nM/M, kon/koff 复合单位, 1/(µM·nm·s)…）。
* **Jacobian 稀疏模式**：构造稀疏依赖，断言传入 `solve_ivp` 的 `jac_sparsity.shape == (n, n)` 且对角为 True。

### 2) A1 对齐门（validate_surrogate）

* 先**忽略**脚本的超时问题（按你的建议），直接在 ad‑hoc `simulate_frozen_model` 路径下对 `A1` 采集：

  * `tumour_volume_l`, `pd1_occupancy`, `tcell_density_per_ul`, `drug_plasma_molar`, `drug_tumor_molar`
* 设定**目标阈值**：`rel_L2 < 1e-3`，`maxRE < 5e-3`；
* 把指标写入 `artifacts/extended_validation/alignment_metrics_extended.csv` 并在 `docs/alignment_tuning_plan.md` 追加“对齐进度表”。
* **等对齐绿灯后**，再小改 `scripts.validate_surrogate`：

  1. 在每个场景循环前后打印 `len(scheduled_list)`、`t=0` 剂量计数、`#triggers`；
  2. 对长时间未进展的循环设置 `max_iters` 安全阀；
  3. 复用 `simulate_frozen_model` 的 schedule 去重逻辑，避免脚本侧重复注入 0 时刻剂量。

---

## 关键改动清单（按文件）

### `src/offline/segment_integrator.py`

* **新增**：`ScheduledDiscontinuity`, `TriggerEventSpec` dataclass
* **变更**：`run_segmented_integration(...)` 统一处理**定时断点 + 状态触发**；
* **必备**：在每个断点 `record(..., phase="pre"/"post")`；断点后 `np.nextafter(t, +∞)`；封装 **边界错误短路**策略。

> 你目前已经支持 ε‑bump / 稀疏 Jacobian / 聚合同刻剂量，此步主要补**状态触发事件** + **record回调**统一 pre/post。

### `src/offline/snapshot.py`

* **新增**：`build_trigger_specs()`（上文代码片段），只为**状态触发**的事件建 `TriggerEventSpec`。
* **完善**：`jacobian_sparsity()`（你已加）；`_convert_parameter_value/_convert_compartment_value` 改为调 `units.py`。
* **可选**：把抗体等药物的 **分子量**（MW）放入 `parameters` 或 `provenance`，供剂量换算。

### `src/offline/units.py`（**新文件**）

* 实现所有单位换算 API；
* 提供**别名适配**（microliter/µL/uL/microlitre、micrometer^3/µm^3、micromolarity/µM 等）；
* 将 `kon_to_L_per_mol_day` 的推导写成可读公式和注释，替换硬编码常数。

### `src/offline/entities.py`

* **不改主字段**；在 `extras` 中新增：`phase`, `discontinuity_type`, `event_index`, `event_name`, `dose_name`, `target`, `time_key`。
* `ScenarioResult.to_frame()` 自动包含这些列（你现有逻辑已支持）。

### `src/offline/simulation.py`

* **主循环瘦身**：只准备 `schedule` + `triggers`，拼装 `apply_discontinuity/reconcile/record` 回调，再调用 integrator。
* `apply_discontinuity`：

  * 对 **dose**：调用 `model.apply_dose()`（内部调用 `units.py` 做 mg→mol→浓度）；
  * 对 **event**：按 `EventEntry.assignments` 执行赋值（已有）；
  * 每次更新后立刻 `reconcile`。
* `record`：把状态写入 `samples` & `extras`（见上）。

### `src/offline/init_conditions.py`（新）

* `solve_init_to_target_volume(...)` + 单元测试。

### `scripts/validate_surrogate.py`（稍后）

* 增加日志与安全阀；使用 `simulate_frozen_model` 的去重/聚合策略。

---

## 代码片段示例（关键处）

**1) 新的 `record` 回调（simulation.py）**

```python
def make_recorder(result_accum):
    def record(t: float, y: np.ndarray, phase: str, event_index: int | None, event_name: str | None,
               *, kind: str | None = None, target: str | None = None):
        key = _tkey(t)
        result_accum["states"][key] = y.copy()
        result_accum["extras"].setdefault("phase", {})[key] = {"cont":0, "pre":1, "post":2}[phase]
        result_accum["extras"].setdefault("discontinuity_type", {})[key] = {
            None:0, "dose":1, "event_immediate":2, "event_delayed":3
        }[kind]
        if event_index is not None:
            result_accum["extras"].setdefault("event_index", {})[key] = event_index
        if event_name is not None:
            result_accum["extras"].setdefault("event_name", {})[key] = event_name
        if target is not None:
            result_accum["extras"].setdefault("target", {})[key] = target
        result_accum["extras"].setdefault("time_key", {})[key] = key
    return record
```

**2) 剂量换算（snapshot.FrozenModel.apply_dose）内调用 `units.py`**

```python
from .units import amount_to_mol

def apply_dose(self, dose: DoseEntry, amount: float, context: Dict[str, float], state: np.ndarray) -> Dict[str, object]:
    entry = self.species_lookup.get(dose.target) or self.species_name_lookup.get(dose.target)
    # 1) 统一成 “mol”
    mw = context.get(f"MW_{dose.name}", None)  # 或从 parameters/provenance 读取
    amount_mol = amount_to_mol(amount, dose.amount_units or "mole", mw)
    # 2) 浓度态则除以 compartment 容量
    if entry and looks_like_concentration(entry):   # 用 interpreted_dimension/units.py 的判断
        V = context.get(entry.compartment, self.compartments.get(entry.compartment))
        if not V or V == 0.0:
            raise RuntimeError(f"Missing volume for {entry.compartment}")
        delta = amount_mol / V
        new_value = context.get(entry.identifier, 0.0) + delta
        self._apply_target_value(entry.identifier, new_value, context, state)
        return {..., "delta_state_value": delta, "delta_amount_mol": amount_mol}
    else:
        new_value = context.get(entry.identifier if entry else dose.target, 0.0) + amount_mol
        self._apply_target_value(entry.identifier if entry else dose.target, new_value, context, state)
        return {..., "delta_state_value": amount_mol, "delta_amount_mol": amount_mol}
```

---

## 交付顺序建议（可并行）

1. **M1**：补齐状态触发事件 + pre/post + 统一 record（最快收敛对齐语义）。
2. **M2**：`units.py` 重构 + `apply_dose`/`_convert_parameter_value` 全量切换（立刻降低 PK/占有率 delta）。
3. **M3**：初始化例程（减少长暂态对齐偏差）。
4. **M4**：体积动态化（输出插件先行）→ 克隆聚合（输出）→ 再逐步进入 ODE。
5. **M5**：测试+A1 门，最后才回到 `validate_surrogate` 的脚本超时点做**最小改动**。

---

## 验收与目标

* A1 数值门：`rel_L2 < 1e-3` 且 `maxRE < 5e-3`（tumour_volume_l、pd1_occupancy、tcell_density_per_ul、PK 通道）。
* 事件回归：

  * 同刻 3 事件（含一个延迟）→ pre/post 顺序与 SimBiology 一致；
  * 边界零步报错→ 不再挂起，产生一次 pre/post 并前进。
* 单位回归：

  * 覆盖 `1/s, 1/min, 1/(M*s), 1/(µM*nm*s)`, `µL/mL/L`, `µM/nM/M`, `mg/µg/mol` 等全部通过；
  * 剂量审计 `delta_state_value`/`delta_amount_mol` 符合维度。

---
- **风险清单（M2）**：
  1. 2D kon 严格模式需要明确的膜厚/绑定长度（默认仍用 legacy 常数 9.8412890625）。
  2. 剂量若缺药物 MW（mg/µg 输入）会抛错；必须在 parameters/config 中显式提供 `MW_<drug>`。
  3. 派生参数 DAG 若存在循环/缺项，加载将失败；需依 audit 脚本 (`scripts/audit_units.py`) 每次检查。
