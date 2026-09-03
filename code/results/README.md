# v7.8 数值结果

Exp2 由 v7.8.0 源码以 4 线程重新运行。Exp1 和 Exp3--7 的算法、算例与输出均未改变，
因此沿用 v7.7.0 的正式结果，不重复计算。

| 主题 | 结果摘要 |
|---|---|
| 主比较 | adaptive/global-reference 均为 13/13 收敛，循环和 4459/28972 |
| 中心扫描 | `m=1,...,128`；最小有效收敛因子为 0.9417637，对应 `m=38` |
| 离线采样参考 | 设计组平均/最大 gap 8.93%/29.65%，验证组 27.34%/54.20% |
| 稳健性与计时 | seed 和 RHS 均全部收敛；平均总时间 523.477/6518.831 ms |
| 停止消融 | adaptive/fixed-step/fixed-residual 收敛数为 6/6、5/6、5/6 |
| 固定物理加密 | adaptive 循环数 122/119/145，有效因子 0.892425/0.890129/0.909090 |
| 三层初试 | adaptive V-cycle 为 101/237 次，有效因子 0.870836/0.943331 |

## 文件

| 实验 | 文件 |
|---|---|
| Exp1 | `experiment1_two_grid_comparison.txt` |
| Exp2 摘要 | `experiment2_step_scan.txt` |
| Exp2 全部 128 个收敛率 | `experiment2_central_step_scan.csv` |
| Exp2 非单调曲线 | `experiment2_central_step_scan.png` |
| Exp3 | `experiment3_oracle_validation.txt` |
| Exp4 | `experiment4_submission_robustness.txt` |
| Exp5 | `experiment5_stopping_ablation.txt` |
| Exp6 | `experiment6_fixed_physical_refinement.txt` |
| Exp7 | `experiment7_multilevel_pilot.txt` |

Exp2 的 CSV 和曲线只包含步数与有效收敛因子，不报告逐步循环数或扫描表格。其他 TXT
保留各自实验所需的循环数、收敛状态、有效收敛因子及辅助指标。
