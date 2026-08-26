# multiscale_two_grids_iteration v4.3

二维高对比扩散问题中，有限 Krylov 能量插值与两网格非单调收敛研究。

- `research_report.md`：v4.3 研究定位、理论版本综合、算法和完整核心结果。
- `theory.tex`：能量路径、粗空间几何、两网格非单调机制与预算化选择理论。
- `VALIDATION_v4.3.md`：严格构建、单元/回归测试、UBSan 和完整实验记录。
- `code/`：C++17 两网格实现、三个整合实验、测试和 TXT 结果。

v4.3 仍以 `expected_rhs_count=R` 作为唯一策略倾向参数。`R=1` 只评估四个候选并使用
16 次 pilot，明显压缩 setup；`R=256` 使用细密候选、160 次 pilot 和至多两个邻域
细化，偏向较少循环数。算法不读取尺度、对比度、拓扑名称或 oracle 信息，并继续只
保留两网格结构。
