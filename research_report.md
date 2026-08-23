# 研究方案与完整阶段报告 v2.9

## 高对比多尺度扩散问题中的局部化能量插值与高效自适应策略

**对应代码：** `two_grids_iteration` v2.9.0  
**统一基线：** $h=1/128$，$H=1/16$，$H/h=8$，对比度 $\kappa=10^4$，随机种子 1  
**两网格框架：** 一次前向 Gauss--Seidel、Galerkin 粗校正、一次后向 Gauss--Seidel

---

## 摘要

本研究考察高对比多尺度扩散方程中，如何在插值稀疏性、构造成本和两网格收敛之间取得平衡。以全局能量极小插值为参照，研究三类近似：限制基函数支撑的空间局部化、使用有限迭代的代数局部化，以及删除小幅值条目的稀疏局部化。

v2.9 延续局部化能量插值主题，并新增通道场有限步 PCG 密集扫描。全部数值结果只保留四类指标：`P density %`、`Setup ms`、`Total ms` 和 `Cycles`。六组实验共 134 个候选，均在 40000 次循环上限内收敛。主要发现如下：

1. 连续随机场中，扩大支撑几乎不减少循环数，低成本几何插值已经足够；
2. 通道和棋盘场明显需要能量修正，但精确全局能量基并不一定带来最小总时间；
3. 全局基存在大量可删除尾部：阈值 $10^{-4}$ 可把密度由约 96%--98% 降至 8%--16%，循环数基本不变；
4. 固定四层支撑中，PCG 容差很快进入性能平台，过度求精只增加 setup；
5. 通道场中，预算式自适应扩展在约 6.24% 密度下，把总时间由 2642 ms 降至 1241 ms；
6. 以全局方程为目标时存在明显有限步窗口：粗扫描中通道场 PCG 32 步为 210 次循环，而精确全局基为 2416 次循环；
7. 密集扫描进一步把当前最佳循环数定位到 PCG 30 步的 202 次，约 28--44 步形成稳定窗口，之后继续逼近精确全局基反而逐步退化；
8. 所有密集扫描候选的粗矩阵求解均通过独立长双精度稠密 Cholesky 交叉检查，因此该非单调现象不能归因于现有粗矩阵直接求解器的数值失效。

这些结果说明，后续方法不应无条件追求更精确或更大的插值，而应按列判断“继续迭代、扩展支撑、剪枝或停止”，并把判断本身的成本纳入优化目标。

---

## 1. 研究问题与范围

### 1.1 模型问题

在 $\Omega=(0,1)^2$ 上考虑

$$
-\nabla\cdot(a(x)\nabla u)=f,
\qquad u|_{\partial\Omega}=0,
$$

其中

$$
0<a_{\min}\le a(x)\le a_{\max},
\qquad \kappa=\frac{a_{\max}}{a_{\min}}.
$$

规则网格五点离散得到线性系统

$$
Au=b.
$$

将未知量分为细点 F 和粗点 C：

$$
A=
\begin{bmatrix}
A_{FF}&A_{FC}\\
A_{CF}&A_{CC}
\end{bmatrix},
\qquad
P(W)=
\begin{bmatrix}
W\\I
\end{bmatrix}.
$$

研究对象是满足 C 点注入的插值矩阵 $P$，粗矩阵取

$$
A_c=P^TAP.
$$

### 1.2 核心问题

本研究集中回答四个问题：

1. 插值基函数需要多大的空间支撑？
2. 固定支撑上的局部系统需要求解到多精确？
3. 哪些小幅值条目可以删除而不显著损害两网格收敛？
4. 能否沿有限步 PCG 路径自动识别最佳停止窗口，并进一步为不同列分配支撑和迭代预算？

### 1.3 当前边界

当前结论来自二维结构网格、固定几何粗点、固定 GS 光滑器、单个 $H/h$、单个对比度和单个随机种子。它们构成方法设计依据，但尚不能直接外推到三维、非结构网格、多层 AMG 或任意随机介质。

---

## 2. 理论基础

### 2.1 离散算子的正定性

### 定理 1：扩散矩阵对称正定

在正系数和齐次 Dirichlet 边界下，离散矩阵 $A$ 为对称正定矩阵。

**证明。** 设网格边权为 $\alpha_{ij}>0$，则对任意内部节点向量 $v$，

$$
v^TAv=
\sum_{(i,j)}\alpha_{ij}(v_i-v_j)^2
+\sum_{i\sim\partial\Omega}\alpha_{i\partial}v_i^2\ge0.
$$

若等号成立，则相邻内部节点取值相同，且与边界相邻的节点为零。网格图连通，因此 $v=0$。故 $A$ 对称正定，其主子矩阵 $A_{FF}$ 也对称正定。证毕。

### 2.2 全局能量极小插值

定义总列能量

$$
J(W)=\frac12\operatorname{tr}\bigl(P(W)^TAP(W)\bigr).
$$

### 定理 2：全局极小点及能量差

$J$ 的唯一极小点为

$$
W_\star=-A_{FF}^{-1}A_{FC},
$$

且对任意 $W=W_\star+E$，

$$
J(W)-J(W_\star)
=\frac12\operatorname{tr}(E^TA_{FF}E).
$$

**证明。** 展开得

$$
J(W)=\frac12\operatorname{tr}(W^TA_{FF}W)
+\operatorname{tr}(W^TA_{FC})
+\frac12\operatorname{tr}(A_{CC}).
$$

其导数为 $A_{FF}W+A_{FC}$。由于 $A_{FF}$ 正定，$J$ 严格凸，驻点唯一。将 $W_\star+E$ 代回，关于 $E$ 的一次项由

$$
A_{FF}W_\star+A_{FC}=0
$$

消失，剩余项即为结论。证毕。

该定理精确描述了“接近全局能量基”的含义，但没有说明它一定使固定 GS 两网格最快。

### 2.3 空间局部化

对第 $j$ 个基函数给定 F 点支撑 $S_j$，令

$$
V_{S_j}=\{v:\operatorname{supp}(v)\subseteq S_j\}.
$$

支撑受限基 $w_{S_j}$ 是全局基 $w_\star$ 在 $V_{S_j}$ 中的能量最佳逼近。

### 定理 3：能量投影与嵌套单调性

对任意 $v\in V_{S_j}$，

$$
\|w_\star-v\|_{A_{FF}}^2
=\|w_\star-w_{S_j}\|_{A_{FF}}^2
+\|w_{S_j}-v\|_{A_{FF}}^2.
$$

若 $S_j\subseteq T_j$，则

$$
\|w_\star-w_{T_j}\|_{A_{FF}}
\le \|w_\star-w_{S_j}\|_{A_{FF}}.
$$

**证明。** 支撑受限极小化的一阶条件为

$$
z^TA_{FF}(w_\star-w_{S_j})=0,
\qquad z\in V_{S_j}.
$$

取 $z=w_{S_j}-v$ 并展开平方得到第一式。又因 $V_{S_j}\subseteq V_{T_j}$，在更大空间中的最小误差不会增大。证毕。

标准 LOD 理论在稳定分解和局部 Poincaré 条件下可得到指数衰减，但当前节点注入空间尚未验证这些条件。因此，本文只使用严格的投影和单调性结论，不把指数衰减当作已证明事实。

### 2.4 代数局部化

固定支撑后，每列需要求解 SPD 系统

$$
A_{SS}x=b_S.
$$

若有限迭代结果为 $x_m$，方程不平衡为 $r_m=b_S-A_{SS}x_m$，则

$$
\|x-x_m\|_{A_{SS}}^2=r_m^TA_{SS}^{-1}r_m.
$$

对角预条件 CG 满足标准上界

$$
\|x-x_m\|_{A_{SS}}
\le
2\left(
\frac{\sqrt{\kappa_S}-1}{\sqrt{\kappa_S}+1}
\right)^m
\|x-x_0\|_{A_{SS}},
$$

其中

$$
\kappa_S=\kappa(D^{-1/2}A_{SS}D^{-1/2}),
\qquad D=\operatorname{diag}(A_{SS}).
$$

这只保证局部系统误差随迭代改善，并不保证两网格循环数或总时间单调改善。

### 引理 4：有限步方法的图传播范围

若初始激励位于节点集合 $G_0$，则对角 Jacobi 或对角预条件 CG 的第 $m$ 步修正满足

$$
\operatorname{supp}(x_m-x_0)\subseteq\mathcal N_m(G_0),
$$

其中 $\mathcal N_m$ 是矩阵图上的 $m$ 跳邻域。

**证明。** 对角缩放不改变非零图，每次乘 $A_{SS}$ 最多将支撑扩展一跳。Jacobi 迭代和 CG Krylov 向量均由有限次矩阵作用的线性组合构成，故结论成立。证毕。

### 2.5 稀疏局部化

对每个插值列 $p_j$ 采用相对阈值

$$
(\mathcal T_\delta p_j)_i=
\begin{cases}
(p_j)_i,& |(p_j)_i|>\delta\max_k|(p_j)_k|,\\
0,&\text{其他},
\end{cases}
$$

并始终保留 C 点注入。设 $\widetilde P=P+D$，则

$$
\widetilde A_c-A_c
=D^TAP+P^TAD+D^TAD,
$$

因此

$$
\|\widetilde A_c-A_c\|_2
\le
2\|A\|_2\|P\|_2\|D\|_2
+\|A\|_2\|D\|_2^2.
$$

该界说明剪枝扰动随 $D$ 连续变化，但单纯幅值阈值不能保证对比度无关，所以必须用两网格实验判断安全区间。

### 2.6 两网格匹配

精确粗校正算子为

$$
\Pi_A(P)=P(P^TAP)^{-1}P^TA.
$$

### 定理 5：粗校正的最佳逼近性质

$\Pi_A(P)$ 是到 $\operatorname{range}(P)$ 的 $A$-正交投影。对任意误差 $e$，

$$
\|e\|_A^2
=\|(I-\Pi_A)e\|_A^2
+\|\Pi_Ae\|_A^2.
$$

**证明。** 直接计算得 $\Pi_A^2=\Pi_A$，且

$$
P^TA(e-\Pi_Ae)=0.
$$

因此剩余误差与粗空间正交，Pythagoras 恒等式成立。证毕。

两网格误差传播为

$$
E_{TG}(P)=S_b(I-\Pi_A(P))S_f.
$$

它同时依赖粗空间与光滑器。全局能量最小化只优化插值列能量，不直接最小化 $\|E_{TG}(P)\|_A$，因此有限步基可能比精确全局基更适合固定 GS。

### 2.7 有限步 PCG、光滑器匹配与非单调性

设一次稳定的驻定迭代写成 $I-M^{-1}A$，其对称化矩阵为

$$
\widetilde M=M(M+M^T-A)^{-1}M^T.
$$

在 $M+M^T-A$ 正定时，固定维数 $n_c$ 的两网格最优粗空间由广义特征问题

$$
Av_i=\lambda_i\widetilde Mv_i,
\qquad
0<\lambda_1\le\cdots\le\lambda_n
$$

的前 $n_c$ 个特征向量张成。

### 定理 6：固定光滑器下的最优粗空间

在上述假设下，所有 $n_c$ 维粗空间中，两网格能量范数收敛因子的最小值由

$$
\mathcal V_{opt}=\operatorname{span}\{v_1,\ldots,v_{n_c}\}
$$

达到。换言之，真正的最优粗空间依赖 $A$ 和光滑器 $M$，而不只依赖插值能量。

**证明要点。** 对称两网格收敛恒等式可把 $\|E_{TG}\|_A^2$ 写成关于 $\widetilde M$-正交投影误差的极大商。对所有 $n_c$ 维子空间取最小，再应用 Courant--Fischer 极小极大原理，即得到由最小的 $n_c$ 个广义特征向量张成的空间。完整形式参见最优插值与 compatible relaxation 理论。证毕。

全局能量极小插值 $P_\star$ 优化的是 $J(P)$。除非

$$
\operatorname{range}(P_\star)=\mathcal V_{opt},
$$

否则它不必使固定 GS 两网格最快。精确 F-relaxation 是二者可能重合的重要特殊情形；本研究使用的点 GS 不属于该特殊情形。

### 命题 7：能量改进不推出两网格改进

若 $\operatorname{range}(P_\star)$ 不是固定光滑器下的最优粗空间，并且最优空间在所选 C 点上的限制块可逆，则存在满足 C 点注入的插值 $\widehat P$ 使

$$
J(\widehat P)>J(P_\star),
\qquad
\|E_{TG}(\widehat P)\|_A<\|E_{TG}(P_\star)\|_A.
$$

**证明。** $P_\star$ 是 $J$ 的唯一极小点，所以任何不同的可行 $\widehat P$ 都有更高能量。取 $\widehat P$ 的值域为最优粗空间并恢复 C 点注入形式；由于 $P_\star$ 的值域非最优，$\widehat P$ 的两网格收敛因子严格更小。证毕。

对角预条件 PCG 从几何初值 $W_0$ 出发时，第 $m$ 步满足

$$
W_m-W_0\in
\mathcal K_m(D^{-1}A_{FF},D^{-1}R_0),
\qquad
R_0=-A_{FC}-A_{FF}W_0.
$$

因此 $W_m$ 既是对能量极小解的近似，也是由次数 $m$ 控制的多项式滤波和图传播过程。随着 $m$ 增大，$J(W_m)$ 单调趋近最小值，但粗空间 $\operatorname{range}(P_m)$ 可以先靠近 $\mathcal V_{opt}$，再转向并不适合点 GS 的 $\operatorname{range}(P_\star)$。这给出了“有限步窗口”的理论解释，但不保证任意问题都出现窗口；窗口是否存在、位于何处，必须由低成本的光滑器相关观测决定。

### 2.8 三类近似的统一关系

设 $P_S$ 为精确支撑基，$P_{S,m}$ 为同一支撑上的有限步近似，则

$$
\|P_\star-P_{S,m}\|_{A,F}^2
=\|P_\star-P_S\|_{A,F}^2
+\|P_S-P_{S,m}\|_{A,F}^2.
$$

加入剪枝后通常只能使用

$$
\|P_\star-\widetilde P_{S,m,\delta}\|_{A,F}
\le E_{space}+E_{algebra}+E_{drop}.
$$

因此，高效方法的关键不是把三项都压到最小，而是在满足两网格需求的前提下，停止收益已经低于成本的动作。

---

## 3. 实验设计

### 3.1 公共配置

| 项目 | 设置 |
|---|---|
| 细网格 | $h=1/128$，$127^2$ 个内部节点 |
| 粗网格 | $H=1/16$，$15^2$ 个内部粗点 |
| 尺度比 | $H/h=8$ |
| 对比度 | $10^4$ |
| 系数场 | continuous、channel、checker |
| 随机种子 | 1 |
| 光滑器 | 一次前向 GS + 一次后向 GS |
| 求解停止条件 | 相对容差 $10^{-6}$ |
| 最大循环 | 40000 |
| 线程数 | 4 |

制造解为

$$
u_\star(x,y)=
\sin(\pi x)\sin(\pi y)
+0.23\sin(3\pi x)\sin(2\pi y)
+0.11\sin(7\pi x)\sin(5\pi y),
$$

右端项取 $b=Au_\star$。

### 3.2 输出指标

| 指标 | 定义 |
|---|---|
| P density % | 插值矩阵非零条目占对应稠密矩阵条目数的百分比 |
| Setup ms | 插值构造、Galerkin 粗矩阵形成和粗矩阵分解的总时间 |
| Total ms | Setup 加一次右端项求解的总时间 |
| Cycles | 达到统一停止条件所需的两网格循环数 |

Field、Method 和 Parameter 仅用于标识实验行，不属于评价指标。

### 3.3 六组实验

| 编号 | 变化因素 | 候选 |
|---|---|---|
| experiment1 | 空间支撑 | 几何、2/3/4 层局部能量基、全局能量基 |
| experiment2 | 代数精度 | 固定四层支撑上的六个 PCG 容差 |
| experiment3 | 稀疏剪枝 | 九个全局基列相对阈值 |
| experiment4 | 支撑选择 | 两层、统一三层、强连接扩展、预算式自适应扩展 |
| experiment5 | 有限步构造 | Jacobi/PCG 的 4/8/16/32/64 步及精确全局参考 |
| experiment6 | 通道密集扫描 | PCG 的 16--64 偶数步及精确全局参考；逐行交叉检查粗求解 |

---

## 4. 完整数值结果

以下六张表与 v2.9.0 代码包中的 CSV 完全一致，共 134 行。时间单位为 ms。experiment1--5 保留 v2.8 基线结果；experiment6 是本版新增的独立完整运行，因此不同表之间的毫秒级时间不用于细小差异比较。

### 4.1 experiment1：支撑半径

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | geometric | P_G | 1.3950 | 4.218 | 43.629 | 87 |
| continuous | local-energy | layers=2 | 6.0922 | 61.634 | 114.611 | 86 |
| continuous | local-energy | layers=3 | 12.5471 | 142.231 | 199.178 | 86 |
| continuous | local-energy | layers=4 | 20.5466 | 271.560 | 346.589 | 86 |
| continuous | global-energy | layers=inf | 95.7222 | 1744.830 | 1979.923 | 86 |
| channel | geometric | P_G | 1.3950 | 1.516 | 6133.290 | 16142 |
| channel | local-energy | layers=2 | 6.0922 | 90.965 | 2767.845 | 4768 |
| channel | local-energy | layers=3 | 12.5471 | 179.230 | 3098.187 | 4262 |
| channel | local-energy | layers=4 | 20.5466 | 384.787 | 3245.251 | 3742 |
| channel | global-energy | layers=inf | 97.8420 | 2306.043 | 8587.365 | 2416 |
| checker | geometric | P_G | 1.3950 | 1.442 | 6154.071 | 16392 |
| checker | local-energy | layers=2 | 6.0922 | 75.971 | 685.850 | 1475 |
| checker | local-energy | layers=3 | 12.5471 | 184.507 | 328.809 | 217 |
| checker | local-energy | layers=4 | 20.5466 | 315.108 | 501.694 | 217 |
| checker | global-energy | layers=inf | 96.9969 | 1886.739 | 2469.773 | 217 |

continuous 中所有能量基均为 86 次循环，扩大支撑没有求解收益。channel 中循环数随支撑增大而下降，但全局基 setup 过高，总时间反而最大。checker 的三层基已经达到 217 次，继续扩大支撑没有收益。

### 4.2 experiment2：固定四层上的 PCG 容差

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | local-energy-4 | tol=1e-02 | 11.1357 | 80.890 | 129.096 | 86 |
| continuous | local-energy-4 | tol=3e-03 | 15.7941 | 92.998 | 143.636 | 86 |
| continuous | local-energy-4 | tol=1e-03 | 18.4505 | 120.350 | 180.682 | 86 |
| continuous | local-energy-4 | tol=3e-04 | 20.0343 | 124.559 | 191.429 | 86 |
| continuous | local-energy-4 | tol=1e-04 | 20.4964 | 124.478 | 186.096 | 86 |
| continuous | local-energy-4 | tol=1e-10 | 20.5466 | 302.569 | 380.085 | 86 |
| channel | local-energy-4 | tol=1e-02 | 14.8391 | 86.067 | 2571.511 | 3855 |
| channel | local-energy-4 | tol=3e-03 | 18.2800 | 115.697 | 2657.862 | 3757 |
| channel | local-energy-4 | tol=1e-03 | 19.6917 | 124.007 | 2743.392 | 3748 |
| channel | local-energy-4 | tol=3e-04 | 20.1398 | 147.525 | 2868.174 | 3745 |
| channel | local-energy-4 | tol=1e-04 | 20.2876 | 147.437 | 2874.532 | 3745 |
| channel | local-energy-4 | tol=1e-10 | 20.5466 | 340.613 | 3281.375 | 3742 |
| checker | local-energy-4 | tol=1e-02 | 14.2278 | 73.766 | 399.724 | 550 |
| checker | local-energy-4 | tol=3e-03 | 17.8336 | 94.448 | 238.840 | 217 |
| checker | local-energy-4 | tol=1e-03 | 19.4081 | 111.316 | 265.056 | 217 |
| checker | local-energy-4 | tol=3e-04 | 19.9827 | 125.084 | 272.762 | 217 |
| checker | local-energy-4 | tol=1e-04 | 20.1788 | 145.665 | 300.953 | 217 |
| checker | local-energy-4 | tol=1e-10 | 20.5466 | 333.019 | 480.579 | 217 |

continuous 的循环数对全部容差都相同。channel 在 $10^{-3}$ 后已接近平台。checker 在 $3\times10^{-3}$ 时已达到严格参考的 217 次循环。$10^{-10}$ 适合作为参考，不适合作为默认构造精度。

### 4.3 experiment3：全局基剪枝

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | global-pruned | drop=0e+00 | 95.7222 | 1627.466 | 1843.339 | 86 |
| continuous | global-pruned | drop=1e-04 | 16.3928 | 1423.633 | 1476.747 | 86 |
| continuous | global-pruned | drop=1e-03 | 9.0549 | 1413.445 | 1453.698 | 86 |
| continuous | global-pruned | drop=3e-03 | 6.2830 | 1409.074 | 1457.077 | 86 |
| continuous | global-pruned | drop=1e-02 | 3.8930 | 1412.339 | 1442.715 | 86 |
| continuous | global-pruned | drop=2e-02 | 2.8038 | 1409.665 | 1440.536 | 86 |
| continuous | global-pruned | drop=3e-02 | 2.2598 | 1405.734 | 1436.996 | 89 |
| continuous | global-pruned | drop=5e-02 | 1.6492 | 1405.327 | 1457.098 | 155 |
| continuous | global-pruned | drop=1e-01 | 0.9962 | 1404.469 | 1570.664 | 468 |
| channel | global-pruned | drop=0e+00 | 97.8420 | 2236.668 | 7894.336 | 2416 |
| channel | global-pruned | drop=1e-04 | 8.4467 | 1952.558 | 3144.370 | 2412 |
| channel | global-pruned | drop=1e-03 | 5.6410 | 1947.095 | 3010.062 | 2584 |
| channel | global-pruned | drop=3e-03 | 4.4629 | 1945.440 | 2953.636 | 2588 |
| channel | global-pruned | drop=1e-02 | 3.1815 | 1947.925 | 2917.461 | 2589 |
| channel | global-pruned | drop=2e-02 | 2.4593 | 1945.965 | 2850.514 | 2594 |
| channel | global-pruned | drop=3e-02 | 2.0604 | 1947.893 | 3046.201 | 3091 |
| channel | global-pruned | drop=5e-02 | 1.5799 | 1948.319 | 4963.369 | 8689 |
| channel | global-pruned | drop=1e-01 | 0.9935 | 1943.756 | 15168.056 | 31953 |
| checker | global-pruned | drop=0e+00 | 96.9969 | 1956.652 | 2457.160 | 217 |
| checker | global-pruned | drop=1e-04 | 7.7556 | 1667.442 | 1762.680 | 218 |
| checker | global-pruned | drop=1e-03 | 5.4953 | 1666.769 | 1753.852 | 218 |
| checker | global-pruned | drop=3e-03 | 4.4253 | 1671.132 | 1750.930 | 218 |
| checker | global-pruned | drop=1e-02 | 3.1966 | 1663.258 | 1818.126 | 426 |
| checker | global-pruned | drop=2e-02 | 2.4822 | 1664.078 | 2188.690 | 1497 |
| checker | global-pruned | drop=3e-02 | 2.0837 | 1661.766 | 2818.570 | 3210 |
| checker | global-pruned | drop=5e-02 | 1.6022 | 1660.591 | 4767.240 | 9044 |
| checker | global-pruned | drop=1e-01 | 1.0059 | 1660.223 | 15580.245 | 33270 |

$\delta=10^{-4}$ 后，三类场的密度分别降至 16.39%、8.45% 和 7.76%，循环数几乎不变。continuous 可继续剪至约 2.80% 而保持 86 次；channel 在 2% 左右仍较稳定；checker 对阈值更敏感，$10^{-2}$ 已从 217 次增至 426 次。由于 setup 仍包含全局基构造，该实验是稀疏性上限分析，不是最终高效算法。

### 4.4 experiment4：支撑选择策略

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | base-local | layers=2 | 6.0922 | 49.183 | 83.813 | 86 |
| continuous | geometric-layer | 2 -> 3 | 12.5471 | 103.457 | 163.579 | 86 |
| continuous | indicator-strong | K=64 one-shot | 6.1082 | 62.285 | 110.827 | 86 |
| continuous | adaptive-budget | R=8, B=128, q=16 | 6.2364 | 177.051 | 227.120 | 86 |
| channel | base-local | layers=2 | 6.0922 | 54.348 | 2641.911 | 4768 |
| channel | geometric-layer | 2 -> 3 | 12.5471 | 124.103 | 3186.308 | 4262 |
| channel | indicator-strong | K=64 one-shot | 6.4855 | 101.633 | 2322.201 | 3960 |
| channel | adaptive-budget | R=8, B=128, q=16 | 6.2390 | 184.846 | 1240.756 | 2582 |
| checker | base-local | layers=2 | 6.0922 | 48.490 | 646.755 | 1475 |
| checker | geometric-layer | 2 -> 3 | 12.5471 | 116.814 | 224.596 | 217 |
| checker | indicator-strong | K=64 one-shot | 6.4838 | 89.523 | 249.296 | 358 |
| checker | adaptive-budget | R=8, B=128, q=16 | 6.2368 | 175.675 | 264.834 | 219 |

continuous 中最简单的两层基最好。checker 中统一三层总时间最低。channel 中 adaptive-budget 以 6.239% 密度将总时间从两层基的 2642 ms 降至 1241 ms，降低约 53.0%。这表明自适应扩展的价值高度依赖系数拓扑，必须包含低成本退出机制。

### 4.5 experiment5：有限步全局目标

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | geometric | P_G | 1.3950 | 1.916 | 29.833 | 87 |
| continuous | Jacobi-global | m=4 | 2.8455 | 16.190 | 38.315 | 69 |
| continuous | Jacobi-global | m=8 | 4.7096 | 35.791 | 58.920 | 65 |
| continuous | Jacobi-global | m=16 | 9.1300 | 120.120 | 169.291 | 63 |
| continuous | Jacobi-global | m=32 | 20.4547 | 429.260 | 479.934 | 64 |
| continuous | Jacobi-global | m=64 | 48.1166 | 1973.266 | 2070.809 | 67 |
| continuous | PCG-global | m=4 | 2.8455 | 45.993 | 65.914 | 61 |
| continuous | PCG-global | m=8 | 4.7096 | 72.634 | 108.292 | 66 |
| continuous | PCG-global | m=16 | 9.1300 | 139.589 | 189.243 | 103 |
| continuous | PCG-global | m=32 | 20.4547 | 274.264 | 335.394 | 86 |
| continuous | PCG-global | m=64 | 48.1166 | 637.909 | 768.592 | 86 |
| continuous | global-exact | tol=1e-10 | 95.7222 | 1648.725 | 1855.459 | 86 |
| channel | geometric | P_G | 1.3950 | 1.954 | 6784.469 | 16142 |
| channel | Jacobi-global | m=4 | 2.8455 | 17.197 | 4629.088 | 11604 |
| channel | Jacobi-global | m=8 | 4.7096 | 32.019 | 3624.501 | 9092 |
| channel | Jacobi-global | m=16 | 9.1300 | 109.241 | 4168.438 | 7471 |
| channel | Jacobi-global | m=32 | 20.4547 | 360.327 | 4537.855 | 5176 |
| channel | Jacobi-global | m=64 | 48.1166 | 1786.436 | 4702.974 | 1885 |
| channel | PCG-global | m=4 | 2.8455 | 49.565 | 1939.103 | 5341 |
| channel | PCG-global | m=8 | 4.7096 | 83.828 | 1697.014 | 3734 |
| channel | PCG-global | m=16 | 9.1300 | 150.085 | 696.926 | 1125 |
| channel | PCG-global | m=32 | 20.4547 | 268.758 | 412.039 | 210 |
| channel | PCG-global | m=64 | 48.1166 | 551.607 | 954.107 | 299 |
| channel | global-exact | tol=1e-10 | 97.8420 | 2284.719 | 8730.914 | 2416 |
| checker | geometric | P_G | 1.3950 | 1.241 | 5321.773 | 16392 |
| checker | Jacobi-global | m=4 | 2.8455 | 14.141 | 4512.721 | 11949 |
| checker | Jacobi-global | m=8 | 4.7096 | 35.979 | 4242.885 | 9402 |
| checker | Jacobi-global | m=16 | 9.1300 | 106.586 | 3826.757 | 7349 |
| checker | Jacobi-global | m=32 | 20.4547 | 357.379 | 4492.106 | 5095 |
| checker | Jacobi-global | m=64 | 48.1166 | 1701.946 | 3532.932 | 1227 |
| checker | PCG-global | m=4 | 2.8455 | 50.434 | 1812.500 | 4556 |
| checker | PCG-global | m=8 | 4.7096 | 88.622 | 1711.251 | 3184 |
| checker | PCG-global | m=16 | 9.1300 | 154.732 | 629.204 | 962 |
| checker | PCG-global | m=32 | 20.4547 | 305.067 | 456.419 | 211 |
| checker | PCG-global | m=64 | 48.1166 | 680.593 | 984.625 | 217 |
| checker | global-exact | tol=1e-10 | 96.9969 | 1831.299 | 2419.839 | 217 |

continuous 中 Jacobi 4 步只需 38 ms，进一步逼近全局基没有意义。channel 中 PCG 32 步为 210 次、412 ms，明显优于精确全局基的 2416 次、8731 ms。checker 中 PCG 32 步为 211 次、456 ms，PCG 64 步和精确基虽然循环数接近，但总时间更高。有限步最佳窗口是后续自适应停止器的直接研究对象。

### 4.6 experiment6：通道场 PCG 密集扫描

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| channel | PCG-global | m=16 | 9.1300 | 204.579 | 1230.413 | 1125 |
| channel | PCG-global | m=18 | 10.3088 | 241.255 | 1321.036 | 1645 |
| channel | PCG-global | m=20 | 11.6236 | 220.814 | 1276.274 | 1045 |
| channel | PCG-global | m=22 | 13.0007 | 290.543 | 909.926 | 663 |
| channel | PCG-global | m=24 | 14.4398 | 289.587 | 739.185 | 521 |
| channel | PCG-global | m=26 | 15.8038 | 342.962 | 546.079 | 311 |
| channel | PCG-global | m=28 | 17.3034 | 315.731 | 539.167 | 250 |
| channel | PCG-global | m=30 | 18.8537 | 279.045 | 462.105 | 202 |
| channel | PCG-global | m=32 | 20.4547 | 318.679 | 493.727 | 210 |
| channel | PCG-global | m=34 | 21.9609 | 312.278 | 536.408 | 228 |
| channel | PCG-global | m=36 | 23.5999 | 335.806 | 531.031 | 214 |
| channel | PCG-global | m=38 | 25.2786 | 378.152 | 583.299 | 214 |
| channel | PCG-global | m=40 | 26.9969 | 383.219 | 598.273 | 215 |
| channel | PCG-global | m=42 | 28.6040 | 399.297 | 638.547 | 215 |
| channel | PCG-global | m=44 | 30.3389 | 378.087 | 597.294 | 216 |
| channel | PCG-global | m=46 | 32.1029 | 415.554 | 660.470 | 218 |
| channel | PCG-global | m=48 | 33.8960 | 428.194 | 649.423 | 216 |
| channel | PCG-global | m=50 | 35.5642 | 496.910 | 773.041 | 217 |
| channel | PCG-global | m=52 | 37.3532 | 476.813 | 744.105 | 218 |
| channel | PCG-global | m=54 | 39.1612 | 534.827 | 826.520 | 222 |
| channel | PCG-global | m=56 | 40.9881 | 565.310 | 903.373 | 233 |
| channel | PCG-global | m=58 | 42.6795 | 584.084 | 964.269 | 256 |
| channel | PCG-global | m=60 | 44.4826 | 559.990 | 997.019 | 280 |
| channel | PCG-global | m=62 | 46.2950 | 646.216 | 1145.690 | 294 |
| channel | PCG-global | m=64 | 48.1166 | 672.917 | 1142.995 | 299 |
| channel | global-exact | tol=1e-10 | 97.8420 | 2515.068 | 9003.419 | 2416 |

密集扫描表明，循环数不是 PCG 步数的单调函数。当前最低点位于 $m=30$，为 202 次；$m=28$--44 构成相对稳定的有效窗口。$m=18$ 的局部反弹说明仅凭相邻检查点的单次下降不能立即停止，实际决策器需要短暂耐心窗口或保留历史最优候选。$m>50$ 后循环数和 setup 同时上升，继续逼近精确全局基已无收益。

experiment6 对每个候选都用独立的长双精度稠密 Cholesky 重新求解同一个 Galerkin 粗系统，并与生产稀疏 Cholesky 比较；全部候选通过预设精度阈值。交叉检查成本不计入表中时间，生产粗求解代码未修改。因此，2416 次的精确全局结果不是由粗矩阵直接求解失败造成的。

---

## 5. 结果综合与当前结论

### 5.1 主要规律

1. **支撑扩展存在场依赖。** 平滑场无需复杂基，通道场需要针对困难区域扩展，棋盘场更适合统一的中等支撑。
2. **局部求解存在精度平台。** 达到足够精度后，继续迭代不会改善循环数。
3. **全局基具有长尾。** 大量条目可删，但直接构造全局基再剪枝的 setup 代价过高。
4. **全局极小不等于两网格最优。** 固定光滑器下，有限步构造可以形成更匹配的粗空间。
5. **有限步存在宽窗口和局部反弹。** 通道场当前有效区间约为 28--44 步，单点贪心停止不足以可靠识别窗口。
6. **粗矩阵求解不是异常来源。** 独立长双精度交叉检查支持“光滑器--粗空间失配”的解释。
7. **单一固定策略不可能普遍最优。** 自适应方法必须能识别简单场并提前退出，也必须能为少数困难列增加预算。

### 5.2 当前最合理的基线

| 场景 | 当前建议基线 | 原因 |
|---|---|---|
| continuous | 几何或低步 Jacobi | setup 极低，循环数已经稳定 |
| channel | PCG 30--32 步或预算式支撑扩展 | 密集扫描最低为 30 步；28--44 步形成有效窗口 |
| checker | 三层局部基或 PCG 32 步 | 达到约 217 次循环且成本可控 |

这些建议只针对当前 $H/h$ 和对比度，不应固化为最终规则。

---

## 6. 高效自适应方法的后续设计

### 6.1 优化目标

对每个候选插值 $P$，目标是

$$
\min_P\quad
T_{setup}(P)+qT_{solve}(P),
$$

其中 $q$ 是预期右端项数量。单右端项时更偏向低 setup；多右端项时允许为更少循环投入更多构造成本。方法同时满足

$$
\operatorname{density}(P)\le d_{max},
\qquad k(P)\le k_{max},
$$

分别限制插值密度和两网格循环数。

### 6.2 第一阶段：整体 PCG 自适应早停

最先实现一个简单、可解释的整体早停器，不立即进入复杂的逐列决策。以几何插值 $P_0$ 为初值，每次增加 $\Delta m=2$ 或 4 步 PCG，并在检查点 $m_k$ 保存候选 $P_{m_k}$。

对少量固定随机向量 $z_\ell$ 先施加一次 GS，使

$$
e_\ell=S_fz_\ell
$$

集中到光滑器难以消除的误差。定义粗空间匹配代理

$$
\widehat\rho(P)^2=
\max_{\ell=1,\ldots,s}
\frac{\|S_b(I-\Pi_A(P))e_\ell\|_A^2}
     {\|e_\ell\|_A^2}.
$$

它不使用 PCG 方程残差来判断“是否足够精确”，而直接估计当前粗空间和固定 GS 的互补程度。结合已测得的一次循环成本，预测

$$
\widehat T(P_m)=T_{setup}(P_m)+
q\,\widehat k(P_m)\widehat t_{cycle}(P_m),
$$

其中

$$
\widehat k(P_m)=
\left\lceil
\frac{\log(10^{-6})}
{\log(\min(1-\varepsilon,\max(\widehat\rho(P_m),\rho_{min})))}
\right\rceil.
$$

其中 $\rho_{min}>0$ 和 $\varepsilon>0$ 只用于避免对数奇异；若 $\widehat\rho\ge1$，则直接判定该候选不具备预测收敛收益。

算法始终保存 $\widehat T$ 最小的历史候选。考虑到 $m=18$ 的局部反弹，停止条件不是“当前一次没有改善”，而是连续 $p=2$ 或 3 个检查点没有刷新历史最优，或者触及密度、步数和 setup 预算。最终返回历史最佳候选，而不是最后一个候选。

```text
输入：A、C/F 划分、几何初值 P0、步长 Δm、耐心 p
best <- P0
连续执行 Δm 步 PCG，并形成检查点 Pm
用少量 GS 过滤探针估计 rho_hat 和总成本 T_hat
若 T_hat 刷新最优，则保存 Pm 并清零 patience
否则 patience 加一
若 patience=p 或达到密度/步数预算，则返回 best
```

### 6.3 第二阶段：按列分配预算

整体早停稳定后，再允许不同列执行四种动作：

1. **Stop：** 保留当前列；
2. **Iterate：** 对该列增加一小段 PCG 步数；
3. **Expand：** 沿强连接图增加少量支撑节点并热启动；
4. **Prune：** 删除低风险尾部条目。

对探针 $e_\ell$，可按列累计粗校正前后的局部能量贡献，标记造成最大剩余误差的少数列。只有这些列获得额外预算；若标记比例过高，则回退到整体 PCG 步进，避免逐列调度成本失控。

### 6.4 局部指示量与收益评分

对第 $j$ 列保留三个辅助量：

$$
\eta_{alg,j}^2\approx r_j^TA_{S_jS_j}^{-1}r_j,
\qquad
\eta_{space,j}^2=
\sum_{i\in\partial S_j}\frac{|(Ap_j)_i|^2}{A_{ii}},
\qquad
\eta_{drop,j}^2=
\sum_{i\in\mathcal D_j}A_{ii}|p_{ij}|^2.
$$

它们分别判断继续迭代、扩展支撑和剪枝风险，但不能单独替代光滑器代理。对动作 $a$ 定义

$$
\operatorname{score}(a)=
q\,\widehat{\Delta k}(a)\widehat t_{cycle}
-\widehat{\Delta T}_{setup}(a),
$$

只执行正收益动作，并受总密度预算限制。所有代理和指示量只用于内部选择；正式结果仍只报告 P density %、Setup ms、Total ms 和 Cycles。

### 6.5 低成本实现约束

1. PCG 状态跨检查点连续复用，不能为每个 $m$ 从零重算；
2. 探针数量先取 2--4，并固定随机种子以降低判断噪声；
3. 只在检查点形成 Galerkin 粗矩阵，检查间隔随 setup 增长可逐步加大；
4. 始终保存历史最佳 $P$，防止越过有限步窗口；
5. 设置硬密度、最大步数、耐心和选择开销预算；
6. 按列阶段只更新被标记列，并复用强连接图、局部编号和已有 PCG 状态；
7. 连续多个配置都选择相近步数时，可把该区间作为下一层或相邻问题的热启动先验，但仍保留回退检查。

### 6.6 理论研究目标

后续理论工作分为四个层次：

1. 由随机子空间或幂迭代理论，给出 $\widehat\rho(P)$ 逼近主导两网格收敛因子的概率误差界；
2. 刻画 PCG 多项式滤波使 $\operatorname{range}(P_m)$ 接近光滑器最优子空间的条件，并解释窗口宽度与谱间隙、对比度和通道连通性的关系；
3. 证明带耐心和历史最优回退的离散停止器，在代理误差有界时返回的预测总成本不超过离散 oracle 的可控倍数；
4. 在固定 $H/h$、有界检查点数和有界每列预算下，证明 setup 与存储随细网格规模近线性增长。

能量误差上界

$$
\|P_\star-P\|_{A,F}
\le C(\eta_{space}+\eta_{alg}+\eta_{drop})
$$

仍用于控制局部化误差，但不再被当作两网格早停的充分条件。如果无法得到对比度无关常数，应明确给出常数对 $\kappa$ 和通道拓扑的依赖。

---

## 7. 高效自适应方法的测试计划

### 7.1 测试矩阵

| 因素 | 取值 |
|---|---|
| 细网格 | $1/64,1/128,1/256$ |
| $H/h$ | 4、8、16 |
| 对比度 | $10^2,10^4,10^6$ |
| 系数场 | continuous、channel、checker、旋转通道、多通道交叉 |
| 随机种子 | 每类至少 10 个 |
| 光滑器 | 点 GS、加权 Jacobi、块 GS、精确 F-relaxation（机制对照） |
| 右端项数量 $q$ | 1、10、100 |
| 线程数 | 1、2、4、8 |

### 7.2 对照方法

1. 几何插值；
2. 固定 2/3/4 层局部能量基；
3. 固定 4/8/16/24/28/30/32/36/44/64 步 Jacobi 和 PCG；
4. 整体 PCG 自适应早停；
5. 按列 PCG 预算和当前 adaptive-budget；
6. 全局基、约束能量极小基及其剪枝结果，仅作为参考；
7. 在全部检查点中直接选择最小 Total ms 的离线 oracle；
8. 可获得时加入 Adaptive AMG、Bootstrap AMG 或 Root-Node AMG 实现作为外部基线。

### 7.3 消融实验

| 消融 | 要回答的问题 |
|---|---|
| 用 PCG 残差替代光滑器代理 | 是否越过两网格最佳窗口 |
| 去掉耐心窗口 | 是否在 $m=18$ 一类局部反弹处过早停止 |
| 去掉历史最佳回退 | 是否返回窗口之后的退化候选 |
| 探针数 1/2/4/8 | 判断可靠性与选择开销如何平衡 |
| 去掉代数或边界指示量 | 逐列阶段是否错误分配迭代和支撑预算 |
| 去掉剪枝动作 | 密度是否持续增长 |
| 去掉收益—成本评分 | 收敛改善是否不足以抵消 setup |
| 去掉 PCG 状态复用 | 从零重算增加多少构造成本 |
| 固定 $q=1$ | 是否错误选择只适合多右端项的昂贵基 |

### 7.4 统计方法

每个配置先预热一次，再独立计时 5 次，报告 Setup ms 和 Total ms 的中位数；不同方法使用相同系数种子做配对比较。P density % 和 Cycles 应为确定性结果。跨种子结果报告中位数和四分位距，并单独列出不收敛配置，而不是用均值掩盖失败。

### 7.5 成功判据

| 目标 | 成功标准 |
|---|---|
| 单右端项效率 | Total ms 中位数不高于最佳固定基的 1.10 倍 |
| 多右端项效率 | $q=10,100$ 时优于最佳固定基，并接近离线 oracle |
| 早停质量 | 所选 PCG 步数的 Total ms 不超过离散 oracle 的 1.15 倍 |
| 窗口识别 | 在存在非单调窗口的配置中不返回窗口之后的明显退化候选 |
| 稀疏性 | 在主要测试矩阵中 P density % 受设定预算控制 |
| 收敛稳健性 | 40000 次上限内的成功率至少 99% |
| 跨尺度稳定性 | 网格细化后 Cycles 不出现系统性倍增 |
| 选择开销 | 自适应判断时间不超过 Setup ms 的 15% |
| 可解释性 | 每次动作都能归因于一个指示量和正收益评分 |

### 7.6 否决与回退条件

出现以下任一情况时，应停止继续复杂化当前自适应策略：

1. 选择开销长期超过节省的求解成本；
2. 跨种子排序不稳定，无法优于简单固定规则；
3. 支撑已接近全局但循环仍高，说明问题在粗点或光滑器；
4. 密度预算频繁耗尽，说明逐列局部扩展不适合该场；
5. 多层递归后粗算子快速变稠密。

对应回退顺序为：固定三层基、统一扩层、调整粗点集合、最后再更换光滑器。这样可以避免把所有失败都错误归因于插值精度。

---

## 8. 后续阶段安排

| 阶段 | 主要任务 | 交付结果 |
|---|---|---|
| A：代理校准 | 校准 GS 探针、检查间隔和总成本预测，与离线 oracle 对比 | 代理误差、耐心参数和阈值范围 |
| B：整体早停 | 实现增量 PCG、GS 探针、耐心窗口和历史最佳回退 | 有限步 PCG 自适应早停器 |
| C：逐列预算 | 实现 Iterate、Expand、Prune 的局部收益模型 | 按列预算的自适应原型 |
| D：规模测试 | 完成多网格、对比度、种子和右端项测试矩阵 | 统计结果与失败类型分析 |
| E：多层推广 | 逐层复用决策器并控制密度 | 多层复杂度和收敛评估 |

阶段 A--C 保持固定粗点和固定 GS，以隔离自适应插值本身的贡献；精确 F-relaxation 只作为机制对照。只有当插值策略稳定后，才进入多层和粗点/光滑器联合设计。

---

## 9. 研究价值、局限与结论

能量最小插值、局部化、Krylov 构造和自适应 AMG 都有成熟基础。本研究不把“使用 PCG”本身作为创新，而是进一步研究经典能量插值目标与固定光滑器之间的失配：PCG 迭代精度持续提高时，两网格性能可以先改善、形成窗口、随后退化。

当前结果最清楚地说明三点：第一，逼近全局能量基的精度与固定 GS 两网格性能并非单调对应；第二，通道场的有效窗口在当前配置下约为 28--44 步，最低循环数出现在 30 步；第三，该现象已通过独立粗矩阵求解交叉检查，不能用粗求解器失效解释。下一阶段的核心任务是用少量 GS 过滤探针自动识别窗口，并在额外 setup 成本超过预期收益之前返回历史最佳候选。

在高完成度下，本研究的定位是经典 AMG 插值问题上的机制发现、目标修正和实质性增量创新：不是追求对所有问题都最复杂的插值，而是得到一个在简单问题上迅速退出、在困难通道上找到有限步窗口、在失效时能够转向支撑、粗点或光滑器设计的高效自适应两网格方法。

---

## 参考文献

1. W. L. Wan, T. F. Chan, B. Smith, “An Energy-Minimizing Interpolation for Robust Multigrid Methods,” *SIAM Journal on Scientific Computing*, 2000. <https://doi.org/10.1137/S1064827598334277>
2. P. S. Vassilevski, “General Constrained Energy Minimization Interpolation Mappings for AMG,” *SIAM Journal on Scientific Computing*, 2010. <https://doi.org/10.1137/080726252>
3. L. N. Olson, J. B. Schroder, R. S. Tuminaro, “A General Interpolation Strategy for Algebraic Multigrid Using Energy Minimization,” *SIAM Journal on Scientific Computing*, 2011. <https://doi.org/10.1137/100803031>
4. A. Målqvist, D. Peterseim, “Localization of Elliptic Multiscale Problems,” *Mathematics of Computation*, 2014. <https://doi.org/10.1090/S0025-5718-2014-02868-8>
5. F. Hellman, A. Målqvist, “Contrast Independent Localization of Multiscale Problems,” *Multiscale Modeling & Simulation*, 2018. <https://doi.org/10.1137/16M1100460>
6. C. Janna, A. Franceschini, J. B. Schroder, L. Olson, “Parallel Energy-Minimization Prolongation for Algebraic Multigrid,” *SIAM Journal on Scientific Computing*, 2023. <https://doi.org/10.1137/22M1513794>
7. J. Brannick, F. Cao, K. Kahl, R. Falgout, X. Hu, “Optimal Interpolation and Compatible Relaxation in Classical Algebraic Multigrid,” 2017. <https://arxiv.org/abs/1703.10240>
8. J. Brannick, S. MacLachlan, J. Schroder, B. Southworth, “The Role of Energy Minimization in Algebraic Multigrid Interpolation,” 2019. <https://arxiv.org/abs/1902.05157>
9. M. Brezina, R. Falgout, S. MacLachlan, T. Manteuffel, S. McCormick, J. Ruge, “Adaptive Algebraic Multigrid,” *SIAM Journal on Scientific Computing*, 2006. <https://doi.org/10.1137/040614402>
10. A. Brandt, J. Brannick, K. Kahl, I. Livshits, “Bootstrap AMG,” *SIAM Journal on Scientific Computing*, 2011. <https://doi.org/10.1137/090752973>
