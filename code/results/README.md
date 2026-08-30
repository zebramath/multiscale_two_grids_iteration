# v5.8 数值结果

四份 TXT 均由 5.8.0 最终源码以完整模式、4 线程生成。相对残量容差为 `1e-6`；
adaptive/reference 的循环上限为 20000，geometric 为 30000。

| 文件 | 内容 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 13 个尺寸、对比度和拓扑问题，含两组 256/16 | adaptive/reference 均 13/13 收敛，循环和 4459/28972；geometric 为 4 converged、9 slow-limit |
| `experiment2_finite_path.txt` | 三条归一化有限 PCG 路径 | 能量超额持续下降，循环数在有限区间达到最小后回升 |
| `experiment3_oracle_validation.txt` | 7 个设计问题和 3 个参数冻结后验证问题的 step--2 离线 oracle | adaptive 的设计组平均/最大 gap 为 8.93%/29.65%，验证组为 27.34%/54.20% |
| `experiment4_submission_robustness.txt` | 五 seed、六 RHS、中心问题五次计时 | seed 循环和 327/866；RHS 循环和 469/1507；adaptive 中心 total 中位数 578.406 ms |

结果逐行报告收敛状态、全程有效收敛因子与末端因子。中心计时以预热后的 Q1、中位数
和 Q3 比较 adaptive 与 reference，且两者均达到正式容差后才生成计时表。留出验证
问题在参数冻结后运行，其结果不参与规则调整。
