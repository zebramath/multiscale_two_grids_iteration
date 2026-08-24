# v3.7 数值结果

六份 TXT 均由 3.7.0 最终源码重新生成。每份文件整合一个公开实验入口下的全部内部研究；本目录不包含 CSV。

| 文件 | 关键结果 |
|---|---|
| `experiment1_localization.txt` | 连续场 local layers=4 达到 global 的 90 次循环；通道自适应支撑在约 6.24% 密度下由 6840 降至 3553 次 |
| `experiment2_pcg_path.txt` | 主通道密集扫描最优为 $m=38$、231 次；新增能量、插值/算子复杂度、粗因子填充和 60 步两网格谱代理 |
| `experiment3_adaptive_oracle.txt` | 主问题自适应选 $m=40$、234 次，距最优 1.30%；新增 3×3 pilot 长度与选择松弛初步消融 |
| `experiment4_robustness.txt` | 18 问题全部收敛；自适应总循环 1900，fixed-$m40$ 为 2431，降低 21.84% |
| `experiment5_diagnostics.txt` | 附加组平均/最大 gap 4.96%/13.79%；压力组为 6.73%/28.10% |
| `experiment6_workload.txt` | 分支、曲环与交叉通道上复用选中层次；前三个问题相对 fixed-$m40$ 的 break-even 为 2、7、3 个 RHS，主问题本次单次计时未回收 |

计时为单次 wall-clock 测量。循环数、密度和 oracle gap 可重复用于算法比较；绝对毫秒数应结合机器与负载解释。
