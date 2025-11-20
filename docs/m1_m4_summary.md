# M1–M4 状态总览（Hybrid 架构结论）

## 当前架构
- **Hybrid**：Python `switches.py` 通过微步长/蓄水池（alignment driver）去 “修补+驱动” MATLAB/SimBiology 导出的 `FrozenModel`。
- 三大白盒（PD‑1、T‐cell、Geometry）已在 Python 端实现，并可在 A1 单剂量场景下跑通；多剂量（A6）仍受限于时间尺度失配。

## 核心痛点
- **Time‑Scale Mismatch**：SimBiology 求解器倾向 0.5 d 级别大步长；Python 端的 PD‑1 结合/PK 输入需要极高时间分辨率，导致多剂量场景（A6）出现超时、回溯伪影、数值不稳定。
- **阶跃输入冲击**：外层大步长把高浓度直接施加到整段 pending 时间，PD‑1 白盒被迫在 10⁻⁶ s 级别处理巨大阶跃，Chunk/Ramp 仍难完全缓解。

## 现状与尝试
- 已加入 **输入斜坡 (ramp)**、分段求解、pending gate/step chunk tunables，并放宽 PD‑1 容差；A1 可靠，A6 仍需耗费大量子步（长时间排空 pending）。
- 最近的 “fast‑pass” 设置（`ramp_chunk_days≈step_chunk_days=0.5`, `rtol=1e-6`, `atol=1e-8`）可启动计算，但在 900 s 内仍未完成 A6；metrics 最新覆盖到 A6，但数值误差偏大。

## 根本瓶颈
- 运行时充斥 dict 查找 + while 循环 + 副作用，无法 JIT / Autodiff，也无法利用 GPU；精度和性能同时受限，阻碍贝叶斯推断和批量仿真。

## 结论与下一步（指向 M5）
- M1–M4 的 Hybrid 修补模式已暴露极限：大步长 → 阶跃输入 → 微步长爆炸 → 超时/不稳定。
- **M5 目标**：重写为纯 Python/JAX 原生 QSP 引擎，抛弃外侧 SimBiology 依赖，获得 JIT + Autodiff 能力，并从设计层面统一时间尺度、数据布局与求解器选择。
