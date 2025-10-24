# TODO – SimBiology 语义对齐计划

> 目标：在 **无需 MATLAB** 的纯 Python 环境中，重放 SimBiology 的执行语义，使数值轨迹与 MATLAB 参考在可接受阈值内一致。

---

## 0. 原则（Principles）
- **以 `getequations` 为唯一真相来源**：`Fluxes:` 与 `ODEs:` 的 RHS 不做二次改写，不再用 `S·v` 组装导数。
- **不做全局手工单位换算**：所有单位/量纲以 `configset.json` 与对象自身单位解释；必要的单位转换（如 µm³ → L）只在明确语义点进行。
- **严格复刻时序**：重复赋值（每步 & 事件后）、并发事件按模型顺序、t=0 初始化顺序（初值/变体/初始赋值/剂量）。

---

## 1. 优先级（P0 → P2）

### P0（阻断项，先做）
- [x] 按 `configset.json` 落实 **TimeUnits/StopTime/SolverOptions/CompileOptions(UnitConversion/DimensionalAnalysis)**。
- [x] **单位/量纲判定**：基于 `species.csv.Units` + `CompileOptions.DefaultSpeciesDimension` 生成 `interpreted_dimension`（amount / conc / amount/vol 等）。
- [x] **V_T 等体积量的单位归一**：若重复赋值来自 µm³ 累加，必须按 UnitConversion 语义转换到 **L**（例如 `V_T_L = V_T_um3 * 1e-15`）。
- [x] **ODE 直译**：`f(t,y)` 中先评估 repeated assignment → 评估 `Fluxes:`（可缓存）→ 逐条评估 `ODEs:`。
- [x] **t=0 快检**：逐条比对 `Fluxes:`/`ODEs:` 数值与 MATLAB 冻结；事件触发布尔判定一致。

### P1（核心语义）
- [x] **重复赋值依赖拓扑**：从表达式 AST 建图，拓扑序求值；事件应用后**再次评估**。
- [x] **事件调度**：多触发同时刻按模型顺序执行；事件后重评 repeated assignment；重启积分（`t_event+ε`）。
- [x] **剂量（Dose）**：将 bolus/infusion 映射到事件/外源速率；t=0 剂量先于积分生效。
- [x] **Species 标志**：`BoundaryCondition/Constant/(NonNegative)` 生效；反应不改 boundary 物种；必要时做非负钳制。

### P2（完善）
- [x] **Rate/Algebraic Rules**：rate rule 进入状态导数；algebraic rule 用 root-finding 或事件式约束。
- [x] **诊断与回归**：导出 `equations_eval_t0.csv`（变量名、Python 值、参考值、相对误差）；轨迹 RMSE 门槛。
- [ ] **对照真值机（可选）**：SBML + RoadRunner 作为差异放大镜（以 `getequations` 结果为准）。

---

## 2. 工件契约（Export Contract）

### 2.1 `configset.json`
- `TimeUnits`, `StopTime`
- `SolverType`, `SolverOptions`（`RelTol`, `AbsTol`, `MaxStep`）
- `CompileOptions`（`UnitConversion`, `DimensionalAnalysis`, `DefaultSpeciesDimension`）

### 2.2 `species.csv`（新增/增强列）
- `name`, `compartment`, `initial_value`, `initial_units`
- `units`（当前单位字符串）
- `boundary_condition` (bool), `constant` (bool), `nonnegative` (bool, optional)
- `interpreted_dimension`（派生：amount / conc / amount/vol / other）

### 2.3 其它
- `rules.csv`: `type(initial|repeated|rate|algebraic)`, `target`, `expr`
- `events.csv`: `index_in_model`, `trigger`, `delay`, `assignments[]`
- `doses.csv`: `type(bolus|infusion)`, `target`, `amount/units`, `time/interval/repeat`, `rate/units`
- `variants.csv`: 名称、赋值列表、是否已应用
- `equations.txt`: 原始 `getequations` 全文

---

## 3. Python 执行语义（frozen_model.py）

### 3.1 ODE 函数顺序
1) 解析配置，构建符号表（t, y, params, compartments）。  
2) **Repeated Assignment**：按依赖图求值至固定点（单步内一次）。  
3) 计算并缓存 **Fluxes**。  
4) 按 `ODEs:` 逐条求 RHS，得到 `dy/dt`。  
5) 执行 **Rate/Algebraic Rules**（如有）。  

### 3.2 事件与剂量
- 以 `solve_ivp` 的事件接口或自定义“零点捕捉”触发；同刻并发按 `index_in_model` 顺序执行全部赋值；  
- 赋值后**再次评估 repeated assignment**；以 `y(t_event+)` 作为新起点，`t ← t_event + ε` 继续积分；  
- 剂量：bolus → 事件式瞬时加量；infusion → 事件包围的区间内叠加外源速率或引入临时 state/rule；  
- t=0 剂量与初始赋值按 SimBiology 初始化时序先后应用。

### 3.3 物种标志与约束
- `BoundaryCondition=true`：反应不改该物种；若有 rate rule，仍由其驱动；  
- `Constant=true`：在仿真期间不可被任何反应/规则改写；  
- `NonNegative=true`：步后投影或事件式截断至 0。

---

## 4. 单位/量纲落地

- 内部以 `configset.TimeUnits` 作为时间基准，不再手工把所有参数改成 “per-day/per-min”。  
- 对体积类量：**只在明确需要的表达式上做单位换算**（例如 `V_T_um3 → L` 用 `1e-15` 因子）。  
- 浓度 vs 数量：若 `units` 缺失则按 `DefaultSpeciesDimension` 解释；否则按单位维度识别（amount、amount/vol…）。  
- 以 `ODEs:` 是否出现 `1/compartment` 作为浓度状态的强信号，但不据此二次改写 RHS。

---

## 5. t=0 快检（必须通过后才允许更新参考轨迹）
- [x] 逐条计算 `Fluxes:`，与冻结参考比对（相对误差 ≤ 1e-9 或给定阈值）。  
- [x] 逐条计算 `ODEs:` RHS，与冻结参考比对。  
- [x] 评估所有 `events.trigger`（布尔），一致性通过。  
- [x] 导出 `equations_eval_t0.csv`：`name, python_value, reference_value, rel_err`。

---

## 6. 验收标准（Definition of Done）
- [ ] Example 1–N：关键观测量轨迹的 RMSE、最大相对误差在阈值内（阈值由模型规模与容差配置决定）。  
- [ ] 全量事件时间戳一致（允许 ±事件容差）。  
- [ ] 非负/边界/常数约束无违例（增加断言与日志）。  
- [ ] 代码层单元测试覆盖：表达式解析、依赖拓扑、并发事件顺序、剂量应用、t=0 快检。  

---

## 7. 日志与诊断（建议）
- `--trace-step`：打印每步触发的事件、更新的 repeated assignment 目标及取值；  
- `--trace-units`：打印单位判定与关键转换因子（如 µm³→L）；  
- `--dump-t0`：输出 t=0 的 Flux/ODE/事件断言文件。

---

## 8. 非目标（Out of Scope，本阶段不做）
- 单位系统的全局自动换算（避免与 `getequations` 叠加造成二次缩放）。  
- 复杂代数规则的高级数值稳定性优化（先保证正确性，再优化）。

---

## 9. 收口冲刺（下一阶段）
- [x] **CI 门槛固化**：工作流跑 `python -m scripts.validate_surrogate --replicates 1 --max-rel-err 1e-6`，t=0 诊断失败直接退出。
- [ ] **参考轨迹版本锁**：在 artefacts 中记录 `equations.txt` / `configset.json` SHA，替换需显式 bump。
- [ ] **事件/剂量回归**：补充并发事件、delay ε 续算、触发方向容差等单元测试。
- [ ] **Algebraic Stress Tests**：引入多根/非线性样例，验证求解策略与失败回退日志。
- [ ] **约束监控**：对 NonNegative/Boundary/Constant 违规计数并输出诊断。
- [ ] **Schema 校验**：为导出工件定义 JSON/Pydantic schema，导入时校验版本。
- [ ] **API/CLI 固化**：统一 `simulate_frozen_model` 参数签名，CLI 增加 `--stop-on-t0-fail`。
- [ ] **跨平台验证**：Mac/Linux/Windows 各跑一次验证误差门槛、事件时间戳一致性。
