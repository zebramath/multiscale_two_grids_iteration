# v3.5 数值结果

六份 TXT 均由 3.5.0 源码在最终验证阶段重新生成。每份文件包含一个整合实验的全部内部研究表格；results 目录不包含 CSV。

| 文件 | 重要结果 |
|---|---|
| `experiment1_localization.txt` | 通道 layers 2/3/4/global 为 6840/6050/5275/3227 cycles；残差预算支撑为 3553 |
| `experiment2_pcg_path.txt` | 主通道步长 2 oracle 为 `m=38, 231 cycles`；`m=40` 为 234 |
| `experiment3_adaptive_oracle.txt` | 主问题自适应 `m=40, 234 cycles`，gap 1.30%；六问题平均 gap 7.81% |
| `experiment4_robustness.txt` | 自适应 18/18 收敛，总 cycles 1858；固定 `m=40` 为 2586 |
| `experiment5_diagnostics.txt` | 两组平均 gap 为 4.77% 和 3.51%，最大为 13.79% 和 28.10% |
| `experiment6_workload.txt` | 相对固定 `m=40` 的 break-even 为 34、35、5、74 RHS |

墙钟时间来自单机单次测量，只适合同一实验内的相对比较。
