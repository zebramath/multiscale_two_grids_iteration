# Numerical results retained in v9.2

本目录中的正式结果来自 v9.1 全量运行。v9.2 只精简项目、合并理论稿导言区并更新版本标识，没有改动数值内核；按要求不重复运行实验，因此原始输出中的 `Version: 9.1.0` 被如实保留。后续用 v9.2 源码复现时会输出 `9.2.0`。

| 前缀 | 内容 |
|---|---|
| `experiment1_` | 13 个主比较问题的完整文本表 |
| `experiment2_` | cross/ring 各 128 点谱路径 CSV、真实谱/RHS 双曲线 PNG、endpoint CSV 与汇总 |
| `experiment3_` | fixed-$H$/fixed-$q$ 16 行尺度、传播、复杂度及谱结果 |
| `experiment4_` | 小规模矩阵导出、gap/一阶方向 CSV、汇总及诊断图 |
| `experiment5_` | seed/RHS 稳健性与五次交替计时 |
| `experiment6_` | 几何与 energy endpoint 基线 |
| `experiment7_` | 三层 V-cycle 试验 |

解释注意：

- `rho_TG` 是对称两网格误差算子的 $A$-能量谱半径估计；`rho_eff` 是指定 RHS 残量历史的有效因子，二者不可互换。
- 谱路径使用 200 步冷启动 Lanczos，`spectral_stage_difference` 比较 100 与 200 步 Ritz 值；双尺度表使用 160 步。
- `correction_support_percent` 是中心粗列实际校正支撑；`reachable_support_percent` 是理论半径 $m-1$ 邻域。浮点支撑阈值为 $10^{-14}$。
- setup/solve 时间只适合在同一次运行、同一硬件和构建配置中横向比较。
- `experiment4_local_diagnostic_matrices.csv` 是稠密理论诊断所需的唯一矩阵输入，保留它可使 gap、主角和一阶方向结果独立复现。
