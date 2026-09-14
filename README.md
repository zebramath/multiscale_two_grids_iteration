# Finite energy-minimization paths for two-grid interpolation (v9.2)

本项目研究一个具体但容易被忽略的机制问题：沿着通向 energy-minimizing/ideal interpolation 的有限 Jacobi--PCG 路径，插值能量严格下降时，固定磨光子下的真实两网格谱半径是否也持续改善？答案是否定的。在所研究的高对比度二维扩散类上，有限路径点可显著优于高精度能量端点，并同时保持更低的 setup、插值密度和粗算子复杂度。

v9.2 不再把经验步数规则包装成“自适应算法”。核心定位是理论机制、真实谱路径与可复现实验；固定 $m=c/h$ 仅作为轻量、可调的经验路径坐标。理论稿已改为单文件自包含形式，不再依赖仓库内的自定义样式包。

## 主要结论

- 对单位注入 $P(W)=[W;I]$，energy endpoint 为 $W^*=-A_{FF}^{-1}A_{FC}$，也是本文坐标下的 ideal interpolation。
- 对固定前向/后向 Gauss--Seidel 与精确 Galerkin 粗解，
  \[
  \rho_{\rm TG}(Z)=\left\|(I+ZZ^T)^{-1/2}(T_F-ZT_*)\right\|_2^2.
  \]
  因而能量距离单调缩小不等价于真实两网格目标单调改善。
- 在简单最大奇异值 gap 且 $T_*v\ne0$ 时，energy endpoint 甚至不是全体单位注入插值中的局部两网格最优点；正文给出显式下降方向与余项控制。
- 第 $m$ 步 PCG 校正的支撑不会越过初始残差在 $A_{FF}$ 图上的半径 $m-1$ 邻域。fixed-$H$/fixed-$q$ 实验的 16 个组合均得到 0 个传播违例。
- 二维规则粗化下
  \[
  \kappa(D_{FF}^{-1/2}A_{FF}D_{FF}^{-1/2})
  \lesssim \chi\left(h^2+\frac{1}{q^2\log(2q)}\right)^{-1}.
  \]
  正文分别证明网格容量尺度和对比度一次幂不能一般性改进，不声称所有受限系数几何都达到统一乘积下界。

## 全量结果摘要

中心问题为 $1/h=128,1/H=16,\chi=10^4$。200 步冷启动 $A$-内积 Lanczos/Ritz 路径给出：

v9.2 未改动数值内核并按要求不重复运行实验；`code/results/` 保留原始 v9.1 全量输出及其版本标识，数值内容与本版源码一致。

| 拓扑 | 最小真实 \(\rho_{\rm TG}\) | 所在 \(m\) | 最小 RHS \(\rho_{\rm eff}\) / \(m\) | energy endpoint \(\rho_{\rm TG}\) |
|---|---:|---:|---:|---:|
| cross-channel | 0.949431677 | 36 | 0.941763675 / 38 | 0.996163213 |
| winding-ring | 0.949386203 | 36 | 0.941012129 / 53 | 0.984367546 |

13 个主比较问题中，固定 $m=\operatorname{round}((1/h)/3)$ 为 13/13 收敛、累计 5,975 个循环；高精度 endpoint 为 13/13 收敛、累计 28,972 个循环。五 seed 和六 RHS 稳健性实验均完整收敛。中心问题五次交替顺序计时中，有限点平均 setup/solve/total 为 0.546/0.325/0.872 s，endpoint 为 2.819/9.078/11.896 s。计时只应在同机同构建下比较。

小规模 $16/4$ 直接诊断得到端点 gap $6.6816\times10^{-3}>0$，理论局部区间内一阶方向符号为 28/28 命中；图公式与直接两网格特征值的最大差为 $4.1\times10^{-14}$。该诊断验证理论对象，不用作停止规则。

## 文件结构

- `theory.tex`：单文件自包含的完整理论稿；所需宏包、数学宏和版式配置均在导言区，技术证明置于附录。
- `research_report.md`：假设、实验设计、定量结论与边界的中文研究报告。
- `code/src/`：稀疏线性代数、扩散离散、PCG 路径、两网格/多层循环与谱诊断。
- `code/experiments/`：7 组正式实验。
- `code/results/`：本次全量运行生成的文本、CSV 与图。
- `code/scripts/run_all_experiments.sh`：quick/full 两种复现入口。

## 复现

需要 C++17 编译器、Python 3、NumPy、SciPy 和 Matplotlib。若有 CMake 则自动使用；否则脚本执行直接构建。

```bash
cd code
./scripts/run_all_experiments.sh quick
./scripts/run_all_experiments.sh full
```

可用 `TGI_THREADS`、`TGI_BUILD_DIR`、`TGI_RESULTS_DIR` 和 `TGI_STEP_TIMEOUT_SECONDS` 调整线程、构建目录、结果目录与单步超时。正式结果使用 4 个线程；随机系数与谱初值均固定 seed。

理论稿使用 XeLaTeX 编译：

```bash
xelatex theory.tex
xelatex theory.tex
xelatex theory.tex
```

除 `ctexart` 文档类外，只需 `geometry`、`fancyhdr`、`mathtools`、`amssymb`、`amsthm`、`bm`、`xcolor`、`hyperref` 和 `cleveref`；仓库内不再需要自定义 `.sty` 文件。

## 文献边界

能量极小插值、Krylov/CG 构造和有限能量最小化迭代已有成熟工作；本项目不把它们重新声明为原创算法。Brannick 等所研究的 optimal interpolation 直接最小化特定两层收敛率，并一般区别于 ideal interpolation；二者在特定 $F$-relaxation/经典 AMG 条件下可等价，但不能无条件移植到本文的全变量对称 Gauss--Seidel 设置。准确引用及 DOI 见 `theory.tex`。

## 声明范围

数值结论针对仓库中定义的二维五点高对比度扩散、规则粗点、单位注入、精确粗解和固定对称 Gauss--Seidel 组合。项目没有宣称固定分数对所有 AMG、离散、粗化、磨光子或系数类普适最优；也没有把 RHS 有效因子当成谱半径。
