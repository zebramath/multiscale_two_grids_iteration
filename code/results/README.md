# v4.1 数值结果

三份 TXT 均由 4.1.0 最终源码完整运行生成。

| 文件 | 覆盖范围 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 11 个单因素尺寸/对比度/拓扑问题 | adaptive 2356 次、exact 11074 次，约 4.70 倍；两个高分离问题达到 13.79 和 12.44 倍 |
| `experiment2_finite_path.txt` | 3 个 cross-channel 有限路径 | 能量持续下降，但两网格循环数在有限优势区间后显著回升 |
| `experiment3_oracle_validation.txt` | 7 个问题、step--2 oracle | 平均 gap 0.62%，最大 gap 2.19% |

所有逐问题数据均保存在对应 TXT 中。计时是单次 wall-clock 辅助量；研究结论主要依据
收敛、循环数、候选步数和插值密度。
