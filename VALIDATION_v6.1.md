# multiscale_two_grids_iteration v6.1 验证记录

验证日期：2026-08-31  
数值源码版本：5.10.0（v6.1 按要求不改代码）  
理论与研究方案版本：6.1  
数值基线来源：v5.10 已验证交付件

## 0. v6.1 变更边界

v6.1 仅增强 `theory.tex`、同步更新 `research_report.md` 与根文档；`code/` 及其
`code/results/` 保持 v5.10 字节不变。因数值实现和正式结果没有变化，本记录继承
v5.10 的构建、Sanitizer、完整实验和数值复验结论，同时新增理论结论、交叉引用、
公式结构以及文件哈希一致性检查。

## 1. 验证环境

- Linux 6.18.35，x86_64；
- AMD EPYC 9V74；当前容器提供 9 个逻辑 CPU；
- GCC 13.3.0；
- C++17；
- CMake 在当前环境中不可用，因此正式验证使用脚本提供的等价直接构建路径；
- XeTeX 3.141592653-2.6-0.999995（TeX Live 2023）；当前环境缺少 `ctexart.cls`。

## 2. 构建与测试

以下检查全部通过：

| 检查 | 配置 | 结果 |
|---|---|---|
| 严格构建 | `-O3 -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic -Werror`，并启用编译器支持的扩展告警 | 通过，无告警 |
| `unit_core` | 4 线程构建路径 | 通过 |
| `regression_v510` | 4 线程构建路径 | 通过 |
| quick 验证 | 4 线程，3 个代表问题 | 通过 |
| 串行 quick 验证 | 1 线程，独立构建目录 | 通过，循环数与 4 线程一致 |
| ASan + UBSan `unit_core` | `-O1 -g -fsanitize=address,undefined` | 通过 |
| ASan + UBSan `regression_v510` | `-O1 -g -fsanitize=address,undefined` | 通过 |

Sanitizer 运行设置为遇错立即停止；由于当前容器不支持可靠的 LeakSanitizer 执行路径，
本项关闭 leak 检查，仅记录 AddressSanitizer 与 UndefinedBehaviorSanitizer 结果。

## 3. 四组完整实验

使用 `TGI_THREADS=4` 执行 `code/scripts/run_validation.sh full`，四组实验全部完成：

| 实验 | 结果 |
|---|---|
| 跨尺寸、对比度与拓扑比较 | adaptive/reference 13/13 收敛，循环和 4459/28972；geometric 4/13 收敛 |
| 有限 PCG 路径 | 三条路径均复现能量单调下降和两网格循环数的有限步非单调优势 |
| 窗口受限 step--2 oracle | 7 个设计问题和 3 个冻结后验证问题全部完成 |
| seed、RHS 与重复计时 | 5 seed、6 RHS 全部完成；两种计时方法均达到正式容差 |

中心 128/16 问题的五次预热后计时为：

| 方法 | setup 中位数 ms | solve 中位数 ms | total 中位数 ms | 循环数 |
|---|---:|---:|---:|---:|
| adaptive | 375.720 | 225.587 | 598.639 | 242 |
| global-reference | 1903.654 | 7252.850 | 9156.505 | 3227 |

global-reference/adaptive 的 total 中位数之比为 15.30。绝对计时依赖运行环境，正式
比较以同一环境、交替顺序和 Q1/中位数/Q3 为准。

## 4. 结果一致性

- 实验 1 和实验 2 除版本号外与 v5.9 正式结果逐行一致；
- 实验 3 的数据表与 v5.9 逐行一致，仅将标题和说明收紧为“窗口受限 oracle”；
- 实验 4 的 seed、RHS、循环数和收敛状态与 v5.9 逐行一致，计时按 v5.10 当前环境重新生成；
- quick 结果写入构建目录下的独立 `quick-results`，已验证不会覆盖正式 full 结果。

## 5. v6.1 理论增强与文稿检查

- PCG 单步勾股账本已提升为任意迭代区间的精确能量尾和；
- 图空间投影已补充闭式正交补、两网格残量范数表达和投影距离的谱/Frobenius 精确式；
- 全局稳定界已补充图坐标 Lipschitz 界，以及能量误差到收敛因子和循环数的显式证书；
- 终点有利方向已由充分条件加强为存在性的充要条件；主导左奇异模态的细坐标自动
  正交，因此唯一附加条件收紧为粗响应矩阵满列秩，并给出可直接代入验证的显式方向；
- 局部展开已导出有利/逆向锥中的线性性能界、二次能量界，以及相干射线 PCG 尾的
  可验证严格反单调条件；
- 最优步数双侧界已补充 $hm_{\rm opt}(h)$ 的显式渐近上下极限；
- `theory.tex` 的 `\label`、`\ref`、`\eqref` 完整且无重复标签；
- 所有 `\cite` 均有对应 `\bibitem`，无重复文献键；
- 在临时副本中仅将缺失的 `ctexart` 文档类替换为标准 `article` 后，XeLaTeX 连续两遍
  通过，未出现未定义控制序列、公式环境错误或悬空交叉引用；该检查只验证 TeX/数学
  结构，不作为中文排版替代；
- 随机矩阵核验通过：图空间投影/正交补恒等式、谱与 Frobenius 投影距离、图坐标
  Lipschitz 界、有利方向显式构造，以及单重主导奇异值情形的一阶有限差分均与公式一致；
- 研究方案、理论稿、代码版本、方法名称、容差、循环上限和核心数值均一致；
- `diff -qr` 与逐文件 SHA-256 均确认 v6.1 的 `code/`（包含 `code/results/`）和 v5.10
  完全一致；
- 当前环境因缺少 `ctexart.cls` 无法完成 XeLaTeX 排版，README 已明确中文 LaTeX 依赖。

## 6. 压缩包复验

最终 ZIP 在独立临时目录解包后复核理论文稿静态结构，并逐文件比较 `code/` 与
`code/results/` 的 v5.10 基线 SHA-256。压缩包不包含构建目录、临时日志或 Git 元数据。
