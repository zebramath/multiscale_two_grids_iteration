# v4.3 数值结果

三份 TXT 均由 4.3.0 最终源码以完整模式、4 线程生成。

| 文件 | 覆盖范围 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 11 个尺寸/对比度/拓扑问题 | fast 2850 次、reuse 2382 次、exact 11074 次；fast 相对 reuse 降低 79.7% setup |
| `experiment2_finite_path.txt` | 3 个 cross-channel 有限路径 | 能量持续下降，但循环数可在有限优势区间后显著回升 |
| `experiment3_oracle_validation.txt` | 7 个问题、step--2 oracle | fast 平均/最大 gap 23.86%/111.11%；reuse 为 1.84%/6.81% |

fast 有意以更粗候选和更短 pilot 换取低 setup，不应解释为 near-oracle 策略。所有逐问题
数据均在对应 TXT 中。计时是单次本机 wall-clock 辅助量；主要可重复证据是收敛性、
循环数、候选步数、插值密度和 oracle gap。
