# v4.7 数值结果

四份 TXT 均由 4.7.0 最终源码以完整模式、4 线程生成，正式两网格循环上限为 12000。

| 文件 | 覆盖范围 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 13 个尺寸/对比度/拓扑问题，含 160/16 与 192/16 | fast/reuse/reference 均 13/13 收敛，循环总和为 4240/3303/24383；geometric 为 3 个收敛、10 个 slow-limit、0 个 diverged |
| `experiment2_finite_path.txt` | 3 个 cross-channel 有限路径 | 能量持续下降，循环数在有限优势区间后回升；截断路径点同时报告末端残量和尾部因子 |
| `experiment3_oracle_validation.txt` | 7 个问题、step--2 oracle | fast 平均/最大 gap 为 23.86%/111.11%，reuse 为 1.84%/6.81% |
| `experiment4_submission_robustness.txt` | 5 个系数 seed、6 个迁移 RHS、5 次重复计时 | 三种能量插值方法的 seed 与 RHS 测试全部收敛；reuse 的 RHS 循环总和为 323，reference 为 1507；报告 Q1/中位数/Q3 |

`converged` 表示达到相对欧氏残量 $10^{-6}$；`slow-limit` 表示达到循环上限时尾部仍
收缩；`diverged` 表示残量非有限或末端持续增长。截断行保留循环数、末端残量和尾部
因子，solve/total 时间留空。主实验的时间汇总采用四种方法共同收敛的三个算例。

全部实验统一使用 `global-reference` 表示相对欧氏残量达到 $10^{-10}$ 的全局能量极小
数值参考。循环数、收敛状态和选择预算构成主要证据；固定问题上的重复计时描述当前
工作负载的 setup--solve 取舍。
