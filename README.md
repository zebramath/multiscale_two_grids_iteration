# multiscale_two_grids_iteration v3.7

二维高对比扩散问题的能量插值、有限步全局 PCG 与两网格收敛研究。

- `research_report.md`：当前项目的完整研究方案、理论主线、算法实现和数值结果。
- `theory.tex`：与 v3.7 诊断量和实验设计同步的理论分析。
- `VALIDATION_v3.7.md`：严格编译、单元/回归测试、UBSan 和六组完整实验记录。
- `code/`：C++17 实现、6 个整合实验、单元/回归测试和完整 TXT 结果。

v3.7 在 v3.6 主算法不变的基础上，增加有限 PCG 路径的能量、层次复杂度、
Cholesky 填充和两网格谱代理诊断，并在 experiment3 中加入初步 pilot 长度与
预测松弛消融。
