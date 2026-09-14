# Numerical results

本目录包含正式实验的完整输出。数据生成配置记录在各结果文件和实验源代码中。

| 前缀 | 内容 |
|---|---|
| `experiment1_` | 13 个主比较问题的完整文本表 |
| `experiment2_` | cross/ring 各 128 点谱路径、双曲线图、endpoint 数据与汇总 |
| `experiment3_` | fixed-$H$/fixed-$q$ 的尺度、传播、复杂度及谱结果 |
| `experiment4_` | 小规模矩阵、gap/一阶方向数据、汇总及诊断图 |
| `experiment5_` | seed/RHS 稳健性与五次交替计时 |
| `experiment6_` | 几何插值与 energy endpoint 基线 |
| `experiment7_` | 三层 V-cycle 试验 |

- `rho_TG` 是对称两网格误差算子的 $A$-能量谱半径；`rho_eff` 是指定 RHS 残量历史的有效因子。
- 谱路径使用 200 步冷启动 Lanczos，`spectral_stage_difference` 比较 100 与 200 步 Ritz 值；双尺度表使用 160 步。
- `correction_support_percent` 是中心粗列实际校正支撑；`reachable_support_percent` 是半径 $m-1$ 邻域。浮点支撑阈值为 $10^{-14}$。
- setup/solve 时间用于同一次运行、同一硬件和构建配置内的横向比较。
- `experiment4_local_diagnostic_matrices.csv` 是 gap、主角和一阶方向诊断的矩阵输入。
