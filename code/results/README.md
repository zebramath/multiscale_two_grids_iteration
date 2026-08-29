# v5.4 数值结果

四份 TXT 均由 5.4.0 最终源码以完整模式、4 线程生成。正式两网格相对残量容差为
`1e-6`，循环上限为 20000。

| 文件 | 内容 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 13 个尺寸、对比度和拓扑问题，含两组 256/16 | fast/reuse/reference 均 13/13 收敛，循环总和 4459/4346/28972；geometric 为 3 converged、10 slow-limit |
| `experiment2_finite_path.txt` | 三条归一化有限 PCG 路径 | 能量超额持续下降，循环数在有限区间达到最小后回升 |
| `experiment3_oracle_validation.txt` | 7 个设计问题和 3 个参数冻结后验证问题的 step--2 离线 oracle | 设计组 fast 平均/最大 gap 8.93%/29.65%，reuse 为 4.07%/11.11%；验证组分别为 27.34%/54.20% 和 19.44%/43.62% |
| `experiment4_submission_robustness.txt` | 五 seed、六 RHS、中心问题五次计时 | seed 循环和 327/278/866；RHS 循环和 469/336/1507；fast 中心 total 中位数 552.731 ms |

每一行同时报告收敛状态、全程有效收敛因子与末端因子。主比较不汇总 wall-clock；时间
只在中心 128/16 问题上以预热后的 Q1、中位数和 Q3 比较 fast、reuse 与 reference。
留出验证问题在参数冻结后运行，其结果未用于再次修改策略。
