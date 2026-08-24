# v4.0 数值结果

四份 TXT 均由 4.0.0 最终源码完整运行生成。实验结构围绕有限 PCG 插值与全局能量极小插值的受控比较重新组织；不再保留独立的复杂度/谱扫描、pilot 参数消融、附加诊断和压力测试入口。

| 文件 | 覆盖范围 | 核心结果 |
|---|---|---|
| `experiment1_two_grid_comparison.txt` | 4 组网格、3 个对比度、6 类拓扑 | 自适应 12/12 收敛，累计 1709 次；global-exact 累计 3654 次；平均密度 37.33% 对 95.23% |
| `experiment2_finite_path.txt` | 3 个代表问题的有限 PCG 路径 | 插值能量单调逼近极小值，而两网格循环数在有限步达到更优区间后可再次增加 |
| `experiment3_oracle_validation.txt` | 8 个问题、step-2 离线 oracle | 平均 gap 8.82%，最大 gap 20.00%，用于评价选择质量 |
| `experiment4_multilevel_comparison.txt` | 与主比较相同的 12 个问题，4–5 层 | 独立 V-cycle：自适应 1717 次、exact 3668 次；作为 PCG 预条件子：313 次、350 次；自适应平均 $C_P$ 为 8.65，exact 为 26.27 |

所有逐问题数据都保存在对应 TXT 中；`research_report.md` 对每个实验展示核心子集和汇总结果。计时为单次 wall-clock 测量，只作为同次运行中的辅助量，理论和算法结论以收敛、迭代数、密度与层次复杂度为主。
