# v4.2 数值结果

三份 TXT 均由 4.2.0 最终源码以完整模式生成。

| 文件 | 覆盖范围 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 11 个尺寸/对比度/拓扑问题 | fast 2458 次、reuse 2382 次、exact 11074 次；几何仅 3/11 收敛 |
| `experiment2_finite_path.txt` | 3 个 cross-channel 有限路径 | 能量持续下降，但循环数在有限优势区间后可显著回升 |
| `experiment3_oracle_validation.txt` | 7 个问题、step--2 oracle | fast 平均/最大 gap 4.95%/25.93%；reuse 为 1.84%/6.81% |

所有逐问题数据均在对应 TXT 中。计时是单次本机 wall-clock 辅助量；主要可重复证据
是收敛性、循环数、候选步数、插值密度和 oracle gap。
