# v4.5 数值结果

四份 TXT 均由 4.5.0 最终源码以完整模式、4 线程生成。

| 文件 | 覆盖范围 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 11 个尺寸/对比度/拓扑问题 | fast 2850 次、reuse 2382 次、reference 11074 次；fast 相对 reuse 降低 79.9% setup |
| `experiment2_finite_path.txt` | 3 个 cross-channel 有限路径 | 能量持续下降，但循环数可在有限优势区间后显著回升 |
| `experiment3_oracle_validation.txt` | 7 个问题、step--2 oracle | 检查低 setup 与低循环两档策略相对离线 oracle 的损失 |
| `experiment4_submission_robustness.txt` | 5 个系数 seed、6 个迁移 RHS、5 次重复计时 | 三种方法 seed 测试均 5/5 收敛；reuse 的 RHS 循环总和为 323，reference 为 1507；报告 Q1/中位数/Q3 |

fast 有意以更粗候选和更短 pilot 换取低 setup，不应解释为 near-oracle 策略。
主实验的 `global-exact` 与第四实验的 `global-reference` 都是容差 $10^{-10}$ 的全局
能量极小数值参考，不是解析精确解。所有逐问题
数据均在对应 TXT 中；循环数、收敛性和选择预算是主要证据，wall-clock 只用于策略成本
方向，重复计时也只代表当前机器与固定算例。
