# multiscale_two_grids_iteration v4.2

二维高对比扩散问题中，有限 Krylov 能量插值与两网格非单调收敛研究。

- `research_report.md`：v4.2 研究定位、算法、完整核心结果和结论边界。
- `theory.tex`：精简后的能量路径、两网格错位与可调选择理论。
- `VALIDATION_v4.2.md`：严格构建、单元/回归测试、UBSan 和完整实验记录。
- `code/`：C++17 两网格实现、三个整合实验、测试和 TXT 结果。

v4.2 以 `expected_rhs_count=R` 统一调节选择成本。`R=1` 使用短 pilot、无细化，
偏向低 setup；大 `R` 使用长 pilot 和至多两个邻域候选，偏向较少循环数。算法不读取
尺度、对比度、拓扑名称或 oracle 信息，并继续只保留两网格结构。
