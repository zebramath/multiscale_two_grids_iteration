# v7.6 完整数值结果

七份 TXT 均由 v7.6.0 最终源码以 `full` 模式、4 线程生成。两网格相对残量容差为
`1e-6`；adaptive/global-reference 的循环上限为 20000，geometric 为 30000。
自适应有限步规则使用网格尺度和矩阵对角尺度比，比例与阈值由设计组校准。

## 四个数值主题

| 主题 | 对应结果 | 核心结果 |
|---|---|---|
| 两网格性能与有限路径 | Exp1--2 | adaptive/global-reference 在主比较中均 13/13 收敛，循环和为 4459/28972；三条路径复现能量下降与两网格循环数的有限步最小值 |
| 规则质量、稳健性与代价 | Exp3--5 | 设计组平均/最大 gap 为 8.93%/29.65%，验证组为 27.34%/54.20%；seed 与 RHS 循环和稳定；total 中位数为 753.910/12152.018 ms；停止消融收敛数为 6/6、5/6、5/6 |
| 固定物理系数场加密 | Exp6 | 三层共享节点失配为 0；adaptive 循环数为 122/119/145 |
| 多层初步验证 | Exp7 | adaptive V-cycle 为 101/237 次，global-reference 为 321/770 次 |

## 原始结果文件

| 实验 | 文件 |
|---|---|
| Exp1 | `experiment1_two_grid_comparison.txt` |
| Exp2 | `experiment2_finite_path.txt` |
| Exp3 | `experiment3_oracle_validation.txt` |
| Exp4 | `experiment4_submission_robustness.txt` |
| Exp5 | `experiment5_stopping_ablation.txt` |
| Exp6 | `experiment6_fixed_physical_refinement.txt` |
| Exp7 | `experiment7_multilevel_pilot.txt` |

结果逐行给出收敛状态、最终残量、全程有效收敛因子、末端因子和插值密度。Exp4 的
墙钟时间报告预热后五次重复的 Q1、中位数和 Q3；跨环境比较以循环数、状态、能量与
稀疏度为主要指标。

在 `code` 目录执行：

```bash
./scripts/run_validation.sh full
./scripts/run_validation.sh supplemental
./scripts/run_validation.sh multilevel
```
