# QSP-IO M5 重构工程白皮书 (Engineering Whitepaper)

这份文档是 QSP-IO 项目 **M5 阶段 (Pure JAX Refactor)** 的终极实施规范。它融合了项目现状、实战派专家的反馈以及架构设计的最佳实践。


**版本**: 2.1 (Final Spec)
**优先级**: P0 (Critical)
**预估工时**: 1-2 Sprint (2周)

-----

## 1\. 核心背景与战略目标 (Executive Summary)

### 1.1 现状 (The Problem)

当前架构 (M1-M4) 是一个 **Hybrid（混合）系统**：Python 脚本 (`switches.py`) 试图通过微步长（micro-stepping）去“修补”和“驱动” MATLAB SimBiology 导出的静态模型 (`FrozenModel`)。

  * **核心痛点**：SimBiology 求解器倾向于迈大步（如 0.5 天），而 Python 端的逻辑（如 PD-1 结合、PK 给药）需要极高的时间分辨率。这种 **“时间尺度失配 (Time-Scale Mismatch)”** 导致了 A6 (多剂量) 场景下的计算超时、回溯伪影（Retroactive Artifacts）和数值不稳定。
  * **根本瓶颈**：当前代码充斥着 `dict` 查找、`while` 循环和副作用，**完全无法被 JIT 编译**，也无法进行自动微分（Autodiff）。这意味着无法利用 GPU 加速，也无法进行高效的贝叶斯推断。

### 1.2 M5 的目标 (The Goal)

**彻底抛弃 MATLAB/SimBiology 依赖**，重写为一个 **纯 Python/JAX 原生** 的 QSP 引擎。

  * **性能**：利用 JAX `vmap` 实现 1000+ 虚拟患者的并行仿真（秒级）。
  * **稳定**：利用 Diffrax 的 `PIDController` 和隐式求解器彻底解决刚性问题。
  * **智能**：利用自动微分（`jax.grad`）实现全局敏感性分析（GSA）和 NUTS 贝叶斯采样。

-----

## 2\. 资产清算 (Assets & Liabilities)

在动手写代码前，必须明确“留什么”和“扔什么”。

### ✅ 核心资产 (保留并移植)

1.  **物理方程 (RHS)**：`src/offline/modules/*_whitebox.py` 里的微分方程逻辑是正确的，需逐行翻译为 JAX 代码。
2.  **参数基准**：`parameters/example1_parameters.json` 是经过验证的参数集。
3.  **验收标准**：`artifacts/validation/A1_reference.csv` 是 M5 输出必须重合的金标准。

### 🗑️ 技术债 (坚决删除)

1.  **`src/offline/snapshot.py`**: 不再读取 MATLAB 快照。模型结构将直接写在 Python 代码里。
2.  **`src/offline/modules/switches.py`**: 这是最大的债。里面的 `ramp_chunk`, `pending_dt`, `history_cap` 都是为了修补 SimBiology 问题而生的 Hack，**在 M5 中一行都不要留**。
3.  **`src/offline/units.py`**: 运行时单位换算效率极低。M5 将在参数加载阶段完成一次性换算。

-----

## 3\. 详细实施方案 (Detailed Implementation)

请在根目录新建 `src/m5/` 包。**严禁直接修改 `src/offline/`**。

### Step 1: 参数系统与预处理 (`src/m5/params.py`)

**核心思想**：Solver 内部**无单位**（Dimensionless）。所有单位换算发生在 JSON 加载的那一瞬间。

**代码规范**：

  * 使用 `NamedTuple` 定义参数结构（这是 JAX Pytree，对 JIT 友好）。
  * 所有字段必须是 `float` 或 `jnp.ndarray`。

<!-- end list -->

```python
import jax.numpy as jnp
from typing import NamedTuple
import json

class QSPParams(NamedTuple):
    # --- 动力学参数 (标准化单位: day, L, mol) ---
    # 必须在这里写清楚单位注释，防止后续混淆
    kon_pd1: float        # 1/(M*day)
    koff_pd1: float       # 1/day
    k_growth: float       # 1/day
    k_kill: float         # 1/(cell*day) - 注意这里是 per cell
    vol_capacity: float   # L
    
    # --- PK 参数 ---
    k_cl: float           # L/day (Clearance)
    v_central: float      # L
    v_peripheral: float   # L
    q_exchange: float     # L/day (Distribution)
    
    # --- 治疗参数 ---
    dose_amt: float       # mole (每次给药量)
    
    # --- 外部耦合参数 (预留给代谢/微环境) ---
    lactate_level: float = 1.0  # 无量纲因子或浓度

def load_params(json_path: str) -> QSPParams:
    with open(json_path) as f:
        raw = json.load(f)
    
    # [CRITICAL] 硬编码单位换算逻辑
    # 必须人工核对每一个参数的原始单位
    MW_NIVO = 1.46e5 # g/mol
    
    return QSPParams(
        # 速率：1/s -> 1/day
        kon_pd1 = float(raw.get('kon_PD1_aPD1', 0)) * 86400.0,
        koff_pd1 = float(raw.get('koff_PD1_aPD1', 0)) * 86400.0,
        
        # 体积：uL -> L (SimBio 很多体积是 uL)
        vol_capacity = float(raw.get('vol_tumor_max', 1.0)) * 1e-6,
        v_central = float(raw.get('V_C', 3.0)), # 假设 raw 已经是 L，需核对
        
        # 质量 -> 摩尔
        dose_amt = float(raw.get('dose_mg', 0)) * 1e-3 / MW_NIVO,
        
        # ... 其他参数映射 ...
        # 确保所有字段都被赋值
    )
```

### Step 2: 状态定义与向量场 (`src/m5/ode.py`)

**核心思想**：

1.  **QSSA (拟稳态近似)**：PD-1 结合反应（秒级）比肿瘤生长（天级）快太多。**绝对不要把 PD-1 结合写成 ODE**，否则求解器步长会缩到 $10^{-12}$。直接写代数方程。
2.  **代数内化**：浓度 `conc` 不是状态，是根据 `amt` 和 `vol` 实时算出来的中间变量。

<!-- end list -->

```python
from typing import NamedTuple
import jax
import jax.numpy as jnp
from .params import QSPParams

# 定义状态结构 (Pytree)
# 相比 Flat Array，这能极大减少索引错误
class QSPState(NamedTuple):
    amt_central: jnp.ndarray    # Drug amount (mol)
    amt_peripheral: jnp.ndarray # Drug amount (mol)
    amt_tumor: jnp.ndarray      # Drug amount (mol)
    tumor_cells: jnp.ndarray    # Number of cells
    t_cells: jnp.ndarray        # Number of cells
    # 注意：Volume 和 PD1_Occ 不在状态里，它们是代数变量

def qsp_vector_field(t, y: QSPState, args: QSPParams):
    # 1. 计算中间代数变量 (Algebraic Equations)
    # -------------------------------------------------------
    # 肿瘤体积 (假设线性关系)
    vol_per_cell = 1e-12 # L/cell (示例值)
    vol_tumor = y.tumor_cells * vol_per_cell
    
    # 实时浓度计算 (解决 SimBiology 滞后问题的关键)
    # [TRICK] 使用 jnp.maximum 防止除以零导致的 NaN
    safe_vol_tumor = jnp.maximum(vol_tumor, 1e-12)
    conc_central = y.amt_central / args.v_central
    conc_tumor = y.amt_tumor / safe_vol_tumor
    
    # [CRITICAL] QSSA for PD-1 Occupancy
    # 假设反应瞬间平衡: Occ = Conc / (Conc + Kd)
    kd = args.koff_pd1 / args.kon_pd1
    pd1_occupancy = conc_tumor / (conc_tumor + kd)
    
    # [Feature] 模块化接口：代谢影响
    # 可以在这里插入 lookup table 函数
    kill_modifier = 1.0 / (1.0 + 0.5 * args.lactate_level)

    # 2. 定义导数 (RHS)
    # -------------------------------------------------------
    # PK: Central <-> Peripheral, Central -> Tumor
    flux_cp = args.q_exchange * (conc_central - y.amt_peripheral/args.v_peripheral)
    flux_elim = args.k_cl * conc_central
    # 简化的肿瘤渗透逻辑 (需根据实际方程补全)
    flux_ct = 0.0 
    
    d_amt_c = -flux_elim - flux_cp - flux_ct
    d_amt_p = flux_cp
    d_amt_t = flux_ct
    
    # Tumor Dynamics
    # 杀伤率受 PD1 Occupancy 和 代谢 Modifier 双重影响
    growth = args.k_growth * y.tumor_cells * (1 - vol_tumor / args.vol_capacity)
    killing = args.k_kill * y.t_cells * (1 - pd1_occupancy) * kill_modifier * y.tumor_cells
    
    d_tumor = growth - killing
    d_tcell = 0.0 # 需补全 T 细胞动力学

    # 返回导数结构体 (Diffrax 会自动处理)
    return QSPState(
        amt_central=d_amt_c,
        amt_peripheral=d_amt_p,
        amt_tumor=d_amt_t,
        tumor_cells=d_tumor,
        t_cells=d_tcell
    )
```

### Step 3: 求解器与给药循环 (`src/m5/solver.py`)

**核心思想**：利用 JAX 的 `lax.scan` 将多剂量给药建模为 **“状态重置 (State Reset)”**。

  * **旧做法**：SimBiology 用 Event 修改状态，Python Driver 只能被动响应，导致刚性起步。
  * **新做法**：把仿真切成 `[t0, t1], [t1, t2]...`。每一段积分结束后，手动给 `y.amt_central` 加上 `dose`，作为下一段的初值。这在数学上是完美的。

<!-- end list -->

```python
import jax
import jax.numpy as jnp
from diffrax import diffeqsolve, ODETerm, Kvaerno5, PIDController, SaveAt
from .ode import qsp_vector_field, QSPState
from .params import QSPParams

@jax.jit
def simulate_patient(params: QSPParams, y0: QSPState, dose_times: jnp.ndarray):
    
    # 定义单段积分逻辑 (从一次给药到下一次给药)
    def integrate_interval(carry_y, t_interval):
        t_start, t_end = t_interval
        
        # 求解器配置
        # Kvaerno5: L-stable 隐式求解器，专治刚性方程 (比 Tsit5 稳)
        # PIDController: 自适应步长，自动处理浓度剧变
        term = ODETerm(qsp_vector_field)
        solver = Kvaerno5()
        controller = PIDController(rtol=1e-6, atol=1e-9)
        
        # 积分
        sol = diffeqsolve(
            term, solver, t0=t_start, t1=t_end, dt0=0.01, y0=carry_y,
            args=params,
            stepsize_controller=controller,
            max_steps=50000, # 给足步数，防止刚性区域报错
            saveat=SaveAt(ts=jnp.linspace(t_start, t_end, 100)) # 统一采样输出
        )
        
        # 获取这一段的终点状态
        y_final = jax.tree_map(lambda x: x[-1], sol.ys)
        
        # 【核心魔法】给药重置 (State Reset)
        # 瞬间给药：直接修改 amt_central，不需要积分过程
        new_amt_c = y_final.amt_central + params.dose_amt
        # 使用 NamedTuple 的 _replace 方法更新状态
        y_next = y_final._replace(amt_central=new_amt_c)
        
        # 返回: (下一个区间的初值, 这一段的轨迹)
        return y_next, sol.ys

    # 构造时间段: [(0, 14), (14, 28), ...]
    intervals = jnp.stack([dose_times[:-1], dose_times[1:]], axis=1)
    
    # 启动循环 (编译为一个 XLA Kernel)
    _, all_trajectories = jax.lax.scan(integrate_interval, y0, intervals)
    
    # 拼接所有时间段的结果
    # 输出: QSPState(amt_central=Array[TotalTime, ...], ...)
    flat_traj = jax.tree_map(lambda x: jnp.concatenate(x, axis=0), all_trajectories)
    return flat_traj
```

### Step 4: 可微初始化 (`src/m5/init.py`)

**痛点解决**：SimBio 用 Event Listener 寻找“肿瘤长到 2cm 开始治疗”的时间点，这不可导。
**方案**：使用 **逆向求根 (Root Finding)**。

```python
from jax.scipy.optimize import bisect

def find_start_time(params, target_vol):
    # 定义目标函数：f(t) = V(t) - V_target
    def volume_error_at_t(t):
        # 这里只跑简单的肿瘤生长 ODE (无药状态)
        # y_t = y0 * exp(k * t) ... (如果有解析解最好，没有就积分)
        return calculated_vol - target_vol
    
    # 在 [0, 1000] 天范围内寻找根
    t_start = bisect(volume_error_at_t, 0.0, 1000.0)
    return t_start # 这个 t_start 对 params 是可微的！
```

### Step 5: 贝叶斯推断接口 (`src/m5/inference.py`)

**核心逻辑**：将 **定值容器** 和 **概率分布** 结合。

```python
import numpyro
import numpyro.distributions as dist
from .params import load_params
from .solver import simulate_patient

def model(observed_data, json_path):
    # 1. 加载基准参数 (全是定值)
    base_params = load_params(json_path)
    
    # 2. 定义先验 (Priors) - 仅针对 GSA 筛选出的 Top 5 参数
    # 这里的变量是 Tracer (随机变量)，不是 float
    k_growth_sample = numpyro.sample("k_growth", dist.LogNormal(jnp.log(base_params.k_growth), 0.5))
    
    # 3. 参数注入 (Override)
    # 用随机变量替换掉定值
    # JAX 的多态性允许 NamedTuple 里既有 float 也有 Tracer
    current_params = base_params._replace(
        k_growth = k_growth_sample
    )
    
    # 4. 运行仿真
    traj = simulate_patient(current_params, ...)
    
    # 5. 似然函数
    numpyro.sample("obs", dist.Normal(traj.tumor_cells, 0.1), obs=observed_data)
```

-----

## 4\. 工程验收清单 (The Definition of Done)

请按照此清单验收：

1.  [ ] **参数对齐**：`src/m5/params.py` 打印出的数值与 `artifacts/unit_audit/*.csv` 完全一致。
2.  [ ] **A1 复现**：运行 M5 单次仿真，`tumour_volume` 与 `A1_reference.csv` 的 RMSE \< 1e-3。
3.  [ ] **A6 稳定性**：运行 M5 多剂量仿真，耗时 \< 1秒（GPU），且无 NaN。
4.  [ ] **梯度检查**：`jax.grad(final_tumor)(params)` 能算出非零梯度。

-----

## 5\. 工程落地与新同学上手提示

### 5.1 环境/依赖基线
- Python ≥ 3.10；推荐使用 `uv` 或 `pip-tools` 锁定依赖。
- 必备包：`jax`, `jaxlib`（CPU/GPU 对应版本）、`diffrax`, `equinox`, `numpyro`（如需 HMC/NUTS）、`pydantic` 或 `dataclasses` 仅用于 CLI 配置（勿污染 JAX 树）。
- 设备：GPU/TPU 优先；CPU 仅用于单元测试。确保安装匹配的 jaxlib CUDA/cuDNN 版本。
- 测试框架：`pytest`; 格式化：`ruff`/`black`; 类型检查：`pyright` 或 `mypy`（确保 NamedTuple/Dataclass 仍是 pytree）。

### 5.2 仿真数据与参考文件
- 参考真值：`artifacts/validation/A1_reference.csv`（金标准）；`metrics.csv` 中的 A1/A6 误差指标可用于 sanity check。
- 现有参数：`parameters/example1_parameters.json`；若缺字段，可在 `artifacts/unit_audit/*.csv` 与 `scripts/scenario_registry.py` 中查找剂量/时间表。
- 剂量时间表：现有 A1–A6 配置在 `scripts/scenario_registry.py`；请写一个 `src/m5/schedules.py` 读取这些或自定义 JSON/CSV。
- 初值：不要再依赖 SimBio snapshot。为 M5 直接定义 IC 加载器（可从 A1_reference 的 t=0 行推断 `tumour_volume`, `tcell_density` 等）。

### 5.3 状态/观测对齐
- 状态向量建议包含：`amt_central/amt_peripheral/amt_tumor`, `tumor_cells`, `t_cells`（或分子/质量单位一致化后以摩尔/个体表示）。所有中间量（体积、浓度、占用率）用代数计算，不入状态。
- 可观测量需匹配 legacy CSV 列名：`tumour_volume_l`、`pd1_occupancy`、`tcell_density_per_ul`；必要时提供别名（`tumor_volume_l`）。
- 单位基准：时间=day，体积=L，药量=mol，细胞=counts。所有转换在 `load_params` 阶段完成。

### 5.4 事件/给药实现提示
- 采用“分段积分 + 状态重置”模式：每段终点直接 `amt_central += dose_amt`；避免在 ODE 内部做 if/else 或手工事件。
- Dosing/采样不必同步：使用 `SaveAt(ts=...)` 统一采样点（与 A1/A6 的 time grid 对齐）。
- 若需 PK 多隔室或非线性清除，将其写入 RHS，不要回退到 pending bucket/ramp hack。

### 5.5 验证流水线
- 最小验收：A1 仿真 RMSE(tumour_volume_l) < 1e-3，对照 `A1_reference.csv`。记录 `metrics_m5.csv`。
- 稳定性：A6 多剂量在 GPU 上 <1 s，无 NaN/Inf；`pd1_occupancy` 波形与 A1/A6 参考同级别（rel_L2 < 0.5 可接受作为首版）。
- 梯度检查：对关键参数（`kon/koff`, `k_growth`, `k_kill`）运行 `jax.grad`，确认无 `None/NaN`。

### 5.6 目录建议（保持与旧版隔离）
```
src/m5/
  params.py       # 参数加载/单位换算
  state.py        # QSPState 定义、辅助结构
  ode.py          # qsp_vector_field（代数变量+RHS）
  solver.py       # 分段积分+给药重置；支持批量 (vmap)
  schedules.py    # A1–A6 等剂量表/采样表
  inference.py    # numpyro 接口（可选）
tests/
  test_m5_a1.py   # A1 对齐测试
  test_m5_grad.py # 梯度健全性测试
```

### 5.7 迁移注意事项
- **不要** 复制 `switches.py` 的 pending/ramp/bucket 逻辑；M5 由单一求解器掌控步长。
- PD‑1 占用率推荐 QSSA/代数解；若坚持 ODE，务必限定步长/选择隐式解法并确认与 A1/A6 对齐。
- 任何新的单位换算或经验系数，务必在注释中写明来源（参数文件/论文/SimBio 版本）。

### 5.8 交付物
- 七个核心文件（见 5.6）+ 单元测试 + 生成的 `metrics_m5.csv` 与对照图（可将 `artifacts/validation/A1_reference.csv` 叠加）。
- README_NEW 中的 “What’s Next (M5)” 链接到本文件，并列出运行示例命令（CPU/GPU）。

### 5.9 MATLAB/SimBiology 参考跑法（对照验证用）
- 仍可用现有脚本生成基准数据，便于和 M5 对齐：
  ```bash
  /Volumes/AlbertSSD/Applications/MATLAB_R2023b.app/bin/matlab -batch \
    "cd('/Volumes/AlbertSSD/Program/new/qspiopy'); addpath(fullfile(pwd,'matlab','scripts')); \
     dev_pd1_training_probe('pd1_train_0004');"
  python scripts/dev_pd1_probe_diff.py artifacts/dev/pd1_training_probe_pd1_train_0004.csv
  ```
- A1/A6 参考：`artifacts/validation/A1_reference.csv`（金标准）及 `scripts/scenario_registry.py` 中的时间/剂量表。必要时可用 `scripts/validate_surrogate` 生成新的参考 CSV 作为 M5 对照。

### 5.10 函数式模块化落地建议
- 按机理拆分纯函数，接口统一为 `f(t, y_sub, params) -> dy_sub`，无副作用，便于单测与 vmap：
  ```python
  def pk_module(t, y_pk, params): ...
  def pd1_module(t, y_pd1, params): return d_pd1, h_pd1  # 可返回占用率供其他模块用
  def tumor_module(t, y_tumor, params): ...
  def tcell_module(t, y_tcell, params): ...
  ```
- 在 `qsp_vector_field` 中：
  1) 按预定义索引切片全局状态；  
  2) 计算跨模块代数量（V_T、浓度、H_PD1）；  
  3) 组装各模块参数并调用模块函数；  
  4) 拼接导数。示例骨架：
  ```python
  def qsp_vector_field(t, y, params):
      y_pk, y_pd1, y_tumor, y_tcell = slice_state(y, params.layout)
      v_tumor = calc_volume(y_tumor, y_tcell, params)
      conc_tumor = y_pk[IDX_TUMOR] / jnp.maximum(v_tumor, 1e-12)
      d_pd1, h_pd1 = pd1_module(t, y_pd1, {**params.pd1, "conc_tumor": conc_tumor, "tumor_cells": y_tumor[0]})
      d_tumor = tumor_module(t, y_tumor, {**params.tumor, "teff": y_tcell[0], "h_pd1": h_pd1})
      d_tcell = tcell_module(t, y_tcell, {**params.tcell, "h_pd1": h_pd1})
      d_pk = pk_module(t, y_pk, {**params.pk, "v_tumor": v_tumor})
      return jnp.concatenate([d_pk, d_pd1, d_tumor, d_tcell])
  ```
- 状态索引建议用常量/Enum 统一管理，防止硬编码魔数。新增模块（如 MDSC）只需增加状态切片与导数拼接，不触及现有模块内部。
