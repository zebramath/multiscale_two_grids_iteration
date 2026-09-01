# v7.4 数值结果（沿用 v7.2 完整复跑）

七份 TXT 均由 v7.2.0 最终源码以完整模式、4 线程重新生成，不继承或混用早期版本结果。
相对残量容差为 `1e-6`；adaptive/global-reference 的循环上限为 20000，geometric 为
30000。v7.3 收紧理论和方法定位；v7.4 进一步统一术语、论证和输出文案。两版均未改变
数值算法或实验参数，并按要求原样保留这七份结果。实验 3 的 sampled oracle 仅指预设
窗口内的离散比较基准；gap 不被解释为参数最优性或 near-oracle 证据。

| 文件 | 内容 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 13 个尺寸、对比度和拓扑问题，含两组 256/16 | adaptive/global-reference 均 13/13 收敛，循环和 4459/28972；geometric 为 4 converged、9 slow-limit |
| `experiment2_finite_path.txt` | 三条归一化有限 PCG 路径 | 能量超额持续下降，循环数在有限区间达到最小后回升 |
| `experiment3_oracle_validation.txt` | 7 个设计问题和 3 个参数冻结后验证问题的 step-2 窗口受限离线 sampled oracle | adaptive 的设计组平均/最大 gap 为 8.93%/29.65%，验证组为 27.34%/54.20% |
| `experiment4_submission_robustness.txt` | 五 seed、六 RHS、中心问题五次计时 | seed 循环和 327/866；RHS 循环和 469/1507；adaptive/global-reference 中心 total 中位数 594.021/8727.991 ms |
| `experiment5_stopping_ablation.txt` | 六个问题的停止策略消融 | adaptive/fixed-step/fixed-residual 收敛 6/6、5/6、5/6；记录循环和 1147/21098/25658 |
| `experiment6_fixed_physical_refinement.txt` | 固定物理系数的 32/4、64/8、128/16 嵌套加密 | 共享节点失配 0；adaptive 循环数 122/119/145 |
| `experiment7_multilevel_pilot.txt` | 64/16/8、128/16/8 三层 V-cycle 与首层精确两网格配对 | adaptive V-cycle 为 101/237 次，global-reference 为 321/770 次；adaptive 的 $C_P$ 为 12.5962/13.4528，global-reference 为 42.4913/44.1727 |

结果逐行报告收敛状态、全程有效收敛因子与末端因子。中心计时以预热后的 Q1、中位数
和 Q3 比较 adaptive 与 global-reference，且两者均达到正式容差后才生成计时表。留出验证
问题在参数冻结后运行，其结果不参与规则调整。

完整重现或只重现实验 5--6、实验 7 可分别执行：

```bash
../scripts/run_validation.sh full
../scripts/run_validation.sh supplemental
../scripts/run_validation.sh multilevel
```
