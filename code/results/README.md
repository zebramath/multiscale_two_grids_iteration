# v7.7 完整数值结果

全部结果由 v7.7.0 最终源码以 `full` 模式、4 线程生成。常规两网格相对残量容差为
`1e-6`，循环上限为 20000；Exp2 的 64 点逐步扫描使用 12000 的循环上限。

| 主题 | 结果摘要 |
|---|---|
| 主比较 | adaptive/global-reference 均为 13/13 收敛，循环和 4459/28972 |
| 中心扫描 | `m=1,...,64`；最小值为 `m=38`、231 次、有效因子 0.9417637 |
| 离线采样参考 | 设计组平均/最大 gap 8.93%/29.65%，验证组 27.34%/54.20% |
| 稳健性与计时 | seed 和 RHS 均全部收敛；平均总时间 523.477/6518.831 ms |
| 停止消融 | adaptive/fixed-step/fixed-residual 收敛数为 6/6、5/6、5/6 |
| 固定物理加密 | adaptive 循环数 122/119/145，有效因子 0.892425/0.890129/0.909090 |
| 三层初试 | adaptive V-cycle 为 101/237 次，有效因子 0.870836/0.943331 |

## 文件

| 实验 | 文件 |
|---|---|
| Exp1 | `experiment1_two_grid_comparison.txt` |
| Exp2 完整表 | `experiment2_step_scan.txt` |
| Exp2 结构化数据 | `experiment2_central_step_scan.csv` |
| Exp2 曲线 | `experiment2_central_step_scan.png` |
| Exp3 | `experiment3_oracle_validation.txt` |
| Exp4 | `experiment4_submission_robustness.txt` |
| Exp5 | `experiment5_stopping_ablation.txt` |
| Exp6 | `experiment6_fixed_physical_refinement.txt` |
| Exp7 | `experiment7_multilevel_pilot.txt` |

TXT 结果给出逐例循环数、收敛状态和有效收敛因子；Exp1、Exp2 与 Exp6 还给出尾部因子。
Exp4 的墙钟时间是双方预热后五次交替顺序测量的算术平均值。
