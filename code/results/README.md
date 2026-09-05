# v8.1 数值结果

全部七组实验均由 v8.1.0 最终源码以 4 线程重新运行。正式求解容差为相对欧氏残量
`1e-6`；`global-reference` 的列相对残量容差为 `1e-10`。

| 主题 | 结果摘要 |
|---|---|
| 主比较 | adaptive/global-reference 均为 13/13 收敛，循环和为 4459/28972 |
| 完整路径扫描 | cross 与 winding-ring 分别在扫描区间的 `m=38`、`m=53` 取得最小 $\rho_{\mathrm{eff}}$；两条能量路径均单调下降 |
| 离线采样参考 | 设计组平均/最大 gap 为 8.93%/29.65%，验证组为 27.34%/54.20% |
| 稳健性与计时 | seed 和 RHS 均全部收敛；warm-start 重复计时总时间为 852.797±237.914 / 11446.680±1430.207 ms |
| 停止消融 | adaptive/fixed-step/fixed-residual 收敛数为 6/6、5/6、5/6 |
| 固定物理加密 | adaptive 循环数为 122/119/145，$\rho_{\mathrm{eff}}$ 为 0.892425/0.890129/0.909090 |
| 三层初试 | adaptive V-cycle 为 101/237 次，$\rho_{\mathrm{eff}}$ 为 0.870836/0.943331 |

## 文件

| 实验 | 文件 |
|---|---|
| Exp1 | `experiment1_two_grid_comparison.txt` |
| Exp2 摘要 | `experiment2_step_scan.txt` |
| Exp2 cross 全部 128 点 | `experiment2_cross_channel_path.csv` |
| Exp2 cross 双面板图 | `experiment2_cross_channel_path.png` |
| Exp2 winding-ring 全部 128 点 | `experiment2_winding_ring_path.csv` |
| Exp2 winding-ring 双面板图 | `experiment2_winding_ring_path.png` |
| Exp3 | `experiment3_oracle_validation.txt` |
| Exp4 | `experiment4_submission_robustness.txt` |
| Exp5 | `experiment5_stopping_ablation.txt` |
| Exp6 | `experiment6_fixed_physical_refinement.txt` |
| Exp7 | `experiment7_multilevel_pilot.txt` |

Exp2 的 CSV 列为步数、$J(W_m)=\tfrac12\operatorname{tr}(P_m^\top A_hP_m)$、归一化
能量差和实际残量历程得到的 $\rho_{\mathrm{eff}}$。这里的 $\rho_{\mathrm{eff}}$ 与理论
能量范数两网格因子 $\rho_{\mathrm{TG}}$ 严格区分。
