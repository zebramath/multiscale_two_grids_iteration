# multiscale_two_grids_iteration v7.7 验证记录

数值源码版本：7.7.0  
理论与研究报告版本：v7.7

## 结构与文稿检查

- `research_report.md` 仅保留研究问题、方法定义、实验设置、数值结果和数值结论；定理、
  命题及其证明集中在 `theory.tex`。
- 主结果按尺度、对比度和拓扑拆为三张表；各实验均同时报告循环数与有效收敛因子。
- Exp2 替换为中心 128/16 cross-channel 问题的 `m=1,...,64` 完整扫描，并生成 TXT、CSV
  和 PNG 曲线。
- 所有数值对照表均只保留 adaptive 与 global-reference；初始插值仍作为有限 PCG
  构造的内部输入。
- Exp4 仍采用双方预热后五次交替顺序测量，输出仅保留 setup、solve、total 和循环数的
  算术平均值。
- 旧 Exp2 结果名、旧验证文件及旧版本引用已删除。

## 源码检查

- 七个实验入口、核心头文件、Python 绘图脚本和 shell 脚本均无普通源码注释。
- 删除了几何对照分支、对应循环上限常量、计时分位数函数及不再使用的输出字段。
- `run_validation.sh` 通过 `sh -n`；绘图脚本通过 Python 语法编译检查。
- 全部 C++ 入口由 GCC 13.3.0 以 C++17、`-O3 -DNDEBUG -Wall -Wextra -Wpedantic
  -Werror` 及可用的扩展警告编译通过。

## 完整运行

执行：

```bash
cd code
TGI_BUILD_DIR=build-v77 TGI_RESULTS_DIR=results \
TGI_STEP_TIMEOUT_SECONDS=3600 ./scripts/run_validation.sh full
```

七个程序和 Exp2 绘图步骤均正常结束，最终结果均写入 `code/results`，七份 TXT 的版本
字段均为 7.7.0。

| 实验 | 验证结果 |
|---|---|
| Exp1 | adaptive/global-reference 均 13/13 收敛，累计循环数 4459/28972 |
| Exp2 | 64 个步数全部完成；最小值为 `m=38`、231 次、有效因子 0.9417637 |
| Exp3 | 10 个问题的 adaptive 与 sampled reference 均收敛 |
| Exp4 | seed 为 5/5、RHS 为 6/6 收敛；平均总时间 523.477/6518.831 ms |
| Exp5 | adaptive/fixed-step/fixed-residual 收敛数为 6/6、5/6、5/6 |
| Exp6 | 三层共享节点失配为 0；两种方法全部收敛 |
| Exp7 | 四个三层 V-cycle 组合全部收敛 |

曲线文件为 1560×1344 RGBA PNG，并已人工检查坐标、图例、极小值标记与非单调趋势。

## 理论稿编译条件

环境中的 XeLaTeX 可执行文件存在，但缺少 `ctexart.cls`，因此无法在本环境生成理论稿
PDF。理论源码保留完整文档结构；正式排版需要安装包含 `ctex` 的 TeX Live 中文组件。
