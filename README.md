# multiscale_two_grids_iteration v4.1

二维高对比扩散问题中有限 Krylov 能量插值与两网格非单调收敛研究。

- `research_report.md`：v4.1 研究定位、算法、完整核心结果和结论边界。
- `theory.tex`：有限 PCG 正交尾、Ritz 滤波、两网格谱判据和渐进式选择理论。
- `VALIDATION_v4.1.md`：严格编译、单元/回归测试、UBSan 和三个完整实验记录。
- `code/`：C++17 两网格实现、三个整合实验、测试和完整 TXT 结果。

v4.1 完整移除了多层 V-cycle 和多层预条件 PCG，将研究对象重新收束为自适应有限
PCG、全局精确能量插值和几何插值的受控两网格比较。新的渐进式候选竞速通过广筛、
确认、锚点保护和局部 step--2 细化，将七例 oracle 平均 gap 降至 0.62%。
