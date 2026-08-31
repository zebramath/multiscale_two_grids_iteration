# multiscale_two_grids_iteration v6.2 验证记录

验证日期：2026-08-31  
数值源码版本：6.2.0  
理论与研究方案版本：6.2  
基线来源：v6.1 理论稿与已验证数值实现

## 0. v6.2 变更边界

v6.2 深化局部理论主线，同时删除 v6.1 中未进入后续论证的行列式、Frobenius 投影距离、
最优空间半正定重述、有利/逆向锥和精确共线射线条件。代码算法保持不变；清理仅涉及
6.2 版本同步、把带旧版号的 `regression_v510` 改为稳定名称 `regression_core`，以及
删除测试中的历史性 “reusable hierarchy” 措辞。四份正式结果由 6.2.0 源码重新生成。

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
| `regression_core` | 4 线程构建路径 | 通过 |
| quick 验证 | 4 线程，3 个代表问题 | 通过 |
| 串行 quick 验证 | 1 线程，独立构建目录 | 通过，循环数与 4 线程一致 |
| ASan + UBSan `unit_core` | `-O1 -g -fsanitize=address,undefined` | 通过 |
| ASan + UBSan `regression_core` | `-O1 -g -fsanitize=address,undefined` | 通过 |

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
| adaptive | 382.207 | 213.774 | 594.095 | 242 |
| global-reference | 1883.706 | 7394.038 | 9277.744 | 3227 |

global-reference/adaptive 的 total 中位数之比为 15.62。绝对计时依赖运行环境，正式
比较以同一环境、交替顺序和 Q1/中位数/Q3 为准。

## 4. 结果一致性

- 实验 1--3 的数值表与 v6.1 基线一致，仅版本标识更新为 6.2.0；
- 实验 4 的 seed、RHS、循环数和收敛状态与 v6.1 一致，计时按 v6.2 当前环境重新生成；
- quick 结果写入构建目录下的独立 `quick-results`，已验证不会覆盖正式 full 结果。

## 5. v6.2 理论深化与文稿检查

- 固定插值能量预算下的局部最优问题已被显式求解：证明尖锐一阶常数
  $[\operatorname{tr}(Y^TS^{-1}Y)^{-1}]^{-1/2}$，并给出唯一单位能量最优方向；
- 证明满秩时最佳改善为 $\alpha_E\sqrt{2\varepsilon}+O(\varepsilon)$，秩亏时所有
  一阶改善消失且最佳改善至多为 $O(\varepsilon)$；
- 严格反单调定理现同时量化单步能量下降、$\sqrt{\rho_{TG}}$ 增长和谱因子增长；
- 精确共线射线条件已由渐近相干路径替代：$Z_m=t_mH+o(t_m)$ 且
  $t_{m+1}/t_m\to q<1$ 即可推出最终严格反单调；
- 删除未被主论证继续使用的辅助结论，保留 PCG 能量账本、图空间投影、全局稳定界、
  局部 Pareto 律、反单调尾部和步数尺度窗口这一条闭合主线；
- `theory.tex` 的 `\label`、`\ref`、`\eqref` 完整且无重复标签；
- 所有 `\cite` 均有对应 `\bibitem`，无重复文献键；
- 在临时副本中仅将缺失的 `ctexart` 文档类替换为标准 `article` 后，XeLaTeX 连续两遍
  通过，未出现未定义控制序列、公式环境错误或悬空交叉引用；该检查只验证 TeX/数学
  结构，不作为中文排版替代；
- 随机矩阵核验通过：图空间投影/正交补、图坐标 Lipschitz、有利方向构造、Pareto
  尖锐上界及等号方向、秩亏一阶上界，以及单重主导奇异值下的能量归一有限差分；
- 研究方案、理论稿、代码版本、方法名称、容差、循环上限和核心数值均一致；
- 当前环境因缺少 `ctexart.cls` 无法完成 XeLaTeX 排版，README 已明确中文 LaTeX 依赖。

## 6. 压缩包复验

最终 ZIP 在独立临时目录解包后重新执行文稿静态检查、严格构建、`unit_core` 与
`regression_core`。压缩包不包含构建目录、临时日志或 Git 元数据。
