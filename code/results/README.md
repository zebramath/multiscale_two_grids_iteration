# v7.5 完整数值结果

七份 TXT 均由 v7.5.0 最终源码以 `full` 模式、4 线程生成。两网格相对残量容差为
`1e-6`；adaptive/global-reference 的循环上限为 20000，geometric 为 30000。
sampled oracle 是预设有限窗口中的离线离散参考。自适应有限步规则使用网格尺度和
矩阵对角尺度比，比例与阈值由设计组校准。

| 文件 | 内容 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 13 个尺寸、对比度和拓扑问题 | adaptive/global-reference 均 13/13 收敛，循环和 4459/28972；geometric 为 4 converged、9 slow-limit |
| `experiment2_finite_path.txt` | 三条归一化有限 PCG 路径 | 能量超额持续下降，循环数在有限区间达到最小后回升 |
| `experiment3_oracle_validation.txt` | 7 个设计问题与 3 个参数冻结后验证问题 | adaptive 的设计组平均/最大 gap 为 8.93%/29.65%，验证组为 27.34%/54.20% |
| `experiment4_submission_robustness.txt` | 五 seed、六 RHS、中心问题五次计时 | seed 循环和 327/866；RHS 循环和 469/1507；adaptive/global-reference 的 total 中位数为 570.872/8901.375 ms |
| `experiment5_stopping_ablation.txt` | 六个问题的停止策略消融 | adaptive/fixed-step/fixed-residual 收敛 6/6、5/6、5/6；记录循环和 1147/21098/25658 |
| `experiment6_fixed_physical_refinement.txt` | 固定物理系数的三层嵌套加密 | 共享节点失配 0；adaptive 循环数 122/119/145 |
| `experiment7_multilevel_pilot.txt` | 两个三层 V-cycle 与首层精确两网格配对 | adaptive V-cycle 为 101/237 次，global-reference 为 321/770 次 |

结果逐行给出收敛状态、全程有效收敛因子和末端因子。实验 4 的墙钟时间报告预热后
五次重复的 Q1、中位数和 Q3；循环数、状态、能量与稀疏度是跨环境比较的主要指标。

在 `code` 目录执行：

```bash
./scripts/run_validation.sh full
./scripts/run_validation.sh supplemental
./scripts/run_validation.sh multilevel
```
