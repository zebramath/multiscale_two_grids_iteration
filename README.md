# Finite energy-minimization paths for two-grid interpolation

本项目研究有限 Jacobi--PCG 能量极小化路径上的两网格机制。在二维高对比度扩散问题中，插值能量沿路径严格下降，而固定对称 Gauss--Seidel 两网格谱目标可在有限步达到更优值。有限路径点同时降低 setup 时间、插值密度和粗算子复杂度。

固定 $m=c/h$ 是轻量、可调的路径坐标。项目的核心由精确谱几何、局部方向分析、有限传播尺度和真实谱路径组成。

## 主要结论

- 对单位注入 $P(W)=[W;I]$，能量端点为 $W^*=-A_{FF}^{-1}A_{FC}$，并对应当前坐标下的 ideal interpolation。
- 对固定前向/后向 Gauss--Seidel 与精确 Galerkin 粗解，
  \[
  \rho_{\rm TG}(Z)=\left\|(I+ZZ^T)^{-1/2}(T_F-ZT_*)\right\|_2^2.
  \]
  该公式直接呈现能量几何与两网格谱几何的差异。
- 在简单最大奇异值 gap 且 $T_*v\ne0$ 时，方向 $H=u_F(T_*v)^T$ 在端点附近严格降低 $\rho_{\rm TG}$，给出 energy/ideal endpoint 的局部可改进性。
- 第 $m$ 步 PCG 校正支撑包含于初始残差在 $A_{FF}$ 图上的半径 $m-1$ 邻域。fixed-$H$/fixed-$q$ 的 16 个组合全部达到该传播界。
- 二维规则粗化满足
  \[
  \kappa(D_{FF}^{-1/2}A_{FF}D_{FF}^{-1/2})
  \lesssim \chi\left(h^2+\frac{1}{q^2\log(2q)}\right)^{-1}.
  \]
  容量构造和高对比度边构造分别达到网格尺度与对比度的一次幂。

## 数值结果

中心问题为 $1/h=128,1/H=16,\chi=10^4$。200 步冷启动 $A$-内积 Lanczos/Ritz 路径给出：

| 拓扑 | 最小真实 $\rho_{\rm TG}$ | 所在 $m$ | 最小 RHS $\rho_{\rm eff}$ / $m$ | endpoint $\rho_{\rm TG}$ |
|---|---:|---:|---:|---:|
| cross-channel | 0.949431677 | 36 | 0.941763675 / 38 | 0.996163213 |
| winding-ring | 0.949386203 | 36 | 0.941012129 / 53 | 0.984367546 |

13 个主比较问题中，$m=\operatorname{round}((1/h)/3)$ 全部收敛，累计 5,975 个循环；高精度 endpoint 累计 28,972 个循环。中心问题五次交替计时中，有限点平均 setup/solve/total 为 0.456/0.244/0.699 s，endpoint 为 2.602/8.179/10.781 s。

小规模 $16/4$ 诊断得到 gap $6.6816\times10^{-3}$，理论局部区间内一阶方向符号为 28/28 命中，图公式与直接两网格特征值的最大差为 $4.1\times10^{-14}$。

## 文件结构

- `theory.tex`：自包含理论稿，宏包、数学宏和版式配置均位于导言区，技术证明置于附录。
- `references.bib`：论文写作所需的完整 BibTeX 文献库。
- `research_report.md`：理论链条、实验设计和定量结论。
- `code/src/`：线性代数、扩散离散、PCG 路径、两网格循环与谱诊断。
- `code/experiments/`：6 组正式实验。
- `code/results/`：完整文本、CSV 与图。
- `code/scripts/run_all_experiments.sh`：quick/full 复现入口。

## 复现

需要 C++17 编译器、Python 3、NumPy、SciPy 和 Matplotlib。

```bash
cd code
./scripts/run_all_experiments.sh quick
./scripts/run_all_experiments.sh full
```

`TGI_THREADS`、`TGI_BUILD_DIR` 和 `TGI_RESULTS_DIR` 分别控制线程、构建目录和结果目录。正式结果使用 4 个线程，随机系数与谱初值使用固定 seed。

理论稿使用 XeLaTeX 编译：

```bash
xelatex theory.tex
bibtex theory
xelatex theory.tex
xelatex theory.tex
```

依赖包括 `ctexart`、`geometry`、`fancyhdr`、`mathtools`、`amssymb`、`amsthm`、`bm`、`xcolor`、`hyperref` 和 `cleveref`。

## 文献关系与实验配置

能量极小插值、Krylov/CG 构造和有限能量极小化迭代构成本文的路径基础。Brannick 等的 optimal interpolation 直接优化特定两层收敛率。本文研究单位注入、规则粗点、精确粗解及全变量对称 Gauss--Seidel 下，energy/ideal endpoint 与真实两网格谱目标之间的机制关系。准确引用及 DOI 见 `references.bib`。
