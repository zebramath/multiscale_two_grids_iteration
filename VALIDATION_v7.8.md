# multiscale_two_grids_iteration v7.8 验证记录

数值源码版本：7.8.0  
理论与研究报告版本：v7.8

## 修改核验

- Exp2 固定中心 128/16 cross-channel、对比度 $10^4$ 的问题，扫描范围由
  `m=1,...,64` 扩展为 `m=1,...,128`。
- Exp2 的 TXT 只保留设置和三个关键有效收敛因子，CSV 只含 `m` 与
  `Effective factor`，研究报告只展示非单调曲线，不再提供逐步表格或循环数。
- 曲线采用单面板，完整绘制 128 个收敛率，并标出最小因子位置与 adaptive 位置。
- 主结果的尺度表只保留 cross-channel 网格序列；中心问题分别出现在尺度、对比度和
  拓扑三个子表中；256/16 winding-ring 作为大规模补充结果单独陈述。
- 版本号、README、研究报告、结果说明和理论稿标题已同步为 v7.8。

## 构建与增量运行

- 七个 C++ 实验入口均由 GCC 13.3.0 以 C++17、`-O3 -DNDEBUG -Wall -Wextra
  -Wpedantic -Werror` 及扩展警告严格编译通过。
- `run_validation.sh` 通过 `sh -n`；Python 绘图脚本通过语法编译检查。
- 按“不变实验不用重跑”的要求，仅重新运行 `experiment2_step_scan` 并重新生成 PNG。
- Exp1 和 Exp3--7 的算法、数据与呈现没有数值变化，正式结果沿用 v7.7.0 输出。

Exp2 增量运行覆盖 $m=1,\ldots,128$。最小有效收敛因子为 0.9417637，位于 $m=38$；
adaptive 在 $m=43$ 的有效收敛因子为 0.9444069；global-reference 为 0.9957271。

曲线文件为 1604×984 RGBA PNG，已人工检查坐标、图例、两个标记和完整非单调趋势。

## 理论稿编译条件

环境中的 XeLaTeX 可执行文件存在，但缺少 `ctexart.cls`，因此无法在本环境生成理论稿
PDF。理论源码保留完整文档结构；正式排版需要安装包含 `ctex` 的 TeX Live 中文组件。
