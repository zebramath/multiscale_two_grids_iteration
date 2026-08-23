# 研究方案与完整阶段报告 v2.11

## 高对比多尺度扩散问题中的局部化能量插值与高效自适应策略

**对应代码：** `two_grids_iteration` v2.11.0  
**统一基线：** $h=1/128$，$H=1/16$，$H/h=8$，对比度 $\kappa=10^4$，随机种子 1  
**两网格框架：** 一次前向 Gauss--Seidel、Galerkin 粗校正、一次后向 Gauss--Seidel

---

## 摘要

本研究考察高对比多尺度扩散方程中，如何在插值稀疏性、构造成本和两网格收敛之间取得平衡。以全局能量极小插值为参照，研究三类近似：限制基函数支撑的空间局部化、使用有限迭代的代数局部化，以及删除小幅值条目的稀疏局部化。

v2.11 将六组实验统一改为常数右端项 $b=\mathbf 1$，不再通过预设解析真解制造右端项；三类高对比系数场下的真解一般没有可用的解析表达式。全部数值结果只保留四类指标：`P density %`、`Setup ms`、`Total ms` 和 `Cycles`。六组实验共 134 个候选，其中 124 个在 40000 次循环内达到相对残差 $10^{-6}$，10 个以 `failed@40000` 明确标记。主要发现如下：

1. 连续随机场相对容易，空间局部化从几何插值的 131 次逐步降至约 90 次，但复杂构造在单次求解总时间上并不占优；
2. 通道场和棋盘场在几何插值下均在 40000 次上限内未收敛，说明常数右端项仍能充分暴露高对比难模态；
3. 全局基存在大量可删除尾部：通道场阈值 $10^{-4}$ 将密度由 97.84% 降到 8.45%，循环数从 3227 变为 3226；
4. 固定四层支撑中，PCG 容差进入明显平台；通道场从 $10^{-2}$ 到 $10^{-10}$ 仅由 5449 次降到 5275 次；
5. 通道场中，预算式自适应支撑在约 6.24% 密度下把两层局部基的 6840 次降到 3553 次；
6. 全局目标有限步 PCG 仍存在强烈非单调性：通道场粗扫描 $m=32$ 为 347 次，而全局精确能量基为 3227 次；
7. 密集扫描进一步把当前最小循环数定位到 $m=38$ 的 231 次，约 $m=36$--50 形成低循环平台，继续逼近全局精确基反而显著退化；
8. experiment6 的全部候选均通过独立长双精度稠密 Cholesky 粗解交叉检查，没有发现能够解释该非单调现象的粗矩阵数值求解问题。

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

当前结论来自二维结构网格、固定几何粗点、固定 GS 光滑器、单个 $H/h$、单个对比度、单个随机种子和固定常数右端项 $b=\mathbf1$。它们构成方法设计依据，但尚不能直接外推到三维、非结构网格、多层 AMG、任意随机介质或任意右端项。特别地，表中的 `Cycles` 是当前右端项下的有限步表现，不等同于最坏模态的渐近收敛因子。

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

这里与 LOD（Localized Orthogonal Decomposition）有直接的离散对应。若定义节点型粗插值 $I_Hv=v_C$，则 $\ker I_H=V_F$；全局能量空间 $\operatorname{range}(P_\star)$ 正是 $V_F$ 的 $A$-正交补。因此 $P_\star$ 可视为节点型框架下的理想多尺度基，而把 $w_\star$ 限制到 $S_j$ 上求解，与 LOD 把全局校正子截断到有限层 patch 上具有相同的“全局正交对象局部化”结构。

两者也有关键差别。经典 LOD 通常借助稳定拟插值、局部 Poincaré 不等式等条件证明校正子的指数衰减；当前简单 C/F 节点注入尚未验证这些条件，因此本文不直接移植对比度无关的指数衰减结论。另一方面，有限步 PCG 还给出一种隐式代数局部化：从局部几何初值出发，第 $m$ 步信息至多传播到矩阵图的 $m$ 跳邻域，可与 LOD 的显式 patch 层数形成对应。更重要的是，LOD 主要关注多尺度空间的解逼近，而本研究还要同时匹配固定 GS 光滑器；这正是全局能量最优与两网格最优可能分离的来源。

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

### 命题 8：固定步 PCG 的能量单调性

对每个粗变量列，PCG 第 $m$ 步在仿射 Krylov 空间 $w_0+\mathcal K_m$ 中极小化对应二次能量。由于

$$
w_0+\mathcal K_m\subseteq w_0+\mathcal K_{m+1},
$$

在没有数值击穿时有

$$
J(P_{m+1})\le J(P_m),
\qquad
J(P_m)\downarrow J(P_\star).
$$

**证明。** CG 的极小化空间随步数嵌套，逐列二次能量不增，求和即得结论。证毕。

该命题只涉及 $J(P_m)$。由于 $E_{TG}(P)$ 同时包含光滑器和由 $P$ 决定的粗空间投影，不能从它推出

$$
\|E_{TG}(P_{m+1})\|_A
\le \|E_{TG}(P_m)\|_A.
$$

命题 7 已严格说明能量排序与固定光滑器下的两网格排序可以相反；当前 PCG 密集扫描进一步给出了沿实际 Krylov 路径发生非单调的数值实例。

### 命题 9：给定右端项循环数与最坏模态收敛因子可以给出不同排序

对称前后光滑时，$E_{TG}$ 在 $A$ 内积下自伴。设其 $A$-标准正交特征向量为 $q_i$，特征值为 $\mu_i$，初始误差为

$$
e_0=\sum_i c_iq_i.
$$

则 $k$ 次循环后

$$
\|e_k\|_A^2
=\|r_k\|_{A^{-1}}^2
=\sum_i c_i^2|\mu_i|^{2k},
$$

而最坏模态渐近因子为

$$
\rho_A(E_{TG})=\max_i|\mu_i|.
$$

**证明。** 由 $A$-自伴性可取 $A$-标准正交特征基，且 $e_k=E_{TG}^ke_0=\sum_i c_i\mu_i^kq_i$。利用正交性和 $r_k=Ae_k$ 立即得到等式。证毕。

因此，即使当前常数右端项没有显著激发某个最慢模态，某个候选仍可能在最坏模态因子稍差时用更少循环达到当前停止条件。反之，可复用预条件器的设计不能只按一个右端项选步数。欧氏残差与上述对偶能量残差满足由 $\lambda_{\min}(A)$、$\lambda_{\max}(A)$ 给出的范数等价，但排序仍可能受模态系数影响。后续自适应方法必须明确优化目标：单次给定右端项采用当前残差探针；多右端项或通用求解器采用随机子空间估计的近最坏模态指标。

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

所有实验统一直接取

$$
b=\mathbf1,
$$

对应连续模型中的常数源项 $f\equiv1$。对于随机连续场、通道场和随机棋盘场，真解一般没有解析表达式，因此本研究不再使用制造解误差评价候选，而统一以 $\|r_k\|_2/\|r_0\|_2\le10^{-6}$ 作为求解停止条件。这样也避免了右端项由同一离散算子反向制造而可能引入的特殊模态偏置。

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

### 3.4 v2.11 的数值问题口径

v2.11 只改变统一数值右端项，不改变插值、Galerkin、PCG、Gauss--Seidel 或粗网格直接解算法。六组实验均从头完整运行，文中的所有数值均直接来自 v2.11.0 的对应 CSV。由于 `Setup ms` 与 `Total ms` 是单次墙钟计时，毫秒级小差异只作参考；主要机制判断优先依据循环数、密度和数量级差异。

---

## 4. 完整数值结果

以下六张表与 v2.11.0 代码包中的 CSV 完全一致，共 134 个候选。时间单位为 ms；六个实验均在常数右端项 $b=\mathbf1$ 下重新完整运行。`failed@40000` 表示达到循环上限仍未满足统一停止条件，并非缺失数据。

### 4.1 experiment1：支撑半径

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | geometric | P_G | 1.3950 | 2.918 | 66.027 | 131 |
| continuous | local-energy | layers=2 | 6.0922 | 60.108 | 115.514 | 104 |
| continuous | local-energy | layers=3 | 12.5471 | 139.666 | 192.608 | 91 |
| continuous | local-energy | layers=4 | 20.5466 | 285.609 | 354.571 | 90 |
| continuous | global-energy | layers=inf | 95.7222 | 1593.209 | 1798.406 | 90 |
| channel | geometric | P_G | 1.3950 | 2.278 | 16219.045 | failed@40000 |
| channel | local-energy | layers=2 | 6.0922 | 85.922 | 3584.638 | 6840 |
| channel | local-energy | layers=3 | 12.5471 | 193.244 | 4062.275 | 6050 |
| channel | local-energy | layers=4 | 20.5466 | 324.920 | 4269.254 | 5275 |
| channel | global-energy | layers=inf | 97.8420 | 2097.170 | 8976.755 | 3227 |
| checker | geometric | P_G | 1.3950 | 1.627 | 14088.300 | failed@40000 |
| checker | local-energy | layers=2 | 6.0922 | 87.602 | 1959.464 | 4144 |
| checker | local-energy | layers=3 | 12.5471 | 182.956 | 336.425 | 243 |
| checker | local-energy | layers=4 | 20.5466 | 350.787 | 533.075 | 229 |
| checker | global-energy | layers=inf | 96.9969 | 1767.620 | 2314.558 | 228 |

### 4.2 experiment2：固定四层上的 PCG 容差

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | local-energy-4 | tol=1e-02 | 11.1357 | 84.160 | 136.501 | 111 |
| continuous | local-energy-4 | tol=3e-03 | 15.7941 | 88.578 | 142.356 | 95 |
| continuous | local-energy-4 | tol=1e-03 | 18.4505 | 133.075 | 191.073 | 91 |
| continuous | local-energy-4 | tol=3e-04 | 20.0343 | 145.929 | 206.800 | 90 |
| continuous | local-energy-4 | tol=1e-04 | 20.4964 | 148.874 | 213.587 | 90 |
| continuous | local-energy-4 | tol=1e-10 | 20.5466 | 289.415 | 346.116 | 90 |
| channel | local-energy-4 | tol=1e-02 | 14.8391 | 100.559 | 3177.548 | 5449 |
| channel | local-energy-4 | tol=3e-03 | 18.2800 | 99.799 | 3364.605 | 5296 |
| channel | local-energy-4 | tol=1e-03 | 19.6917 | 117.920 | 3595.122 | 5283 |
| channel | local-energy-4 | tol=3e-04 | 20.1398 | 141.754 | 3662.380 | 5278 |
| channel | local-energy-4 | tol=1e-04 | 20.2876 | 152.816 | 3741.426 | 5277 |
| channel | local-energy-4 | tol=1e-10 | 20.5466 | 329.279 | 3989.239 | 5275 |
| checker | local-energy-4 | tol=1e-02 | 14.2278 | 78.080 | 719.112 | 1186 |
| checker | local-energy-4 | tol=3e-03 | 17.8336 | 95.384 | 259.308 | 265 |
| checker | local-energy-4 | tol=1e-03 | 19.4081 | 113.542 | 263.329 | 231 |
| checker | local-energy-4 | tol=3e-04 | 19.9827 | 126.201 | 276.556 | 229 |
| checker | local-energy-4 | tol=1e-04 | 20.1788 | 142.166 | 307.816 | 229 |
| checker | local-energy-4 | tol=1e-10 | 20.5466 | 319.186 | 467.406 | 229 |

### 4.3 experiment3：全局基剪枝

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | global-pruned | drop=0e+00 | 95.7222 | 1745.703 | 1945.865 | 90 |
| continuous | global-pruned | drop=1e-04 | 16.3928 | 1487.092 | 1547.665 | 90 |
| continuous | global-pruned | drop=1e-03 | 9.0549 | 1478.718 | 1523.397 | 92 |
| continuous | global-pruned | drop=3e-03 | 6.2830 | 1480.349 | 1523.972 | 94 |
| continuous | global-pruned | drop=1e-02 | 3.8930 | 1469.801 | 1513.910 | 100 |
| continuous | global-pruned | drop=2e-02 | 2.8038 | 1475.611 | 1516.805 | 111 |
| continuous | global-pruned | drop=3e-02 | 2.2598 | 1468.609 | 1512.401 | 118 |
| continuous | global-pruned | drop=5e-02 | 1.6492 | 1473.685 | 1540.734 | 206 |
| continuous | global-pruned | drop=1e-01 | 0.9962 | 1468.007 | 1716.977 | 633 |
| channel | global-pruned | drop=0e+00 | 97.8420 | 2179.299 | 9652.533 | 3227 |
| channel | global-pruned | drop=1e-04 | 8.4467 | 1948.641 | 3462.674 | 3226 |
| channel | global-pruned | drop=1e-03 | 5.6410 | 1948.042 | 3430.613 | 3481 |
| channel | global-pruned | drop=3e-03 | 4.4629 | 1946.125 | 3274.254 | 3486 |
| channel | global-pruned | drop=1e-02 | 3.1815 | 1944.195 | 3145.762 | 3488 |
| channel | global-pruned | drop=2e-02 | 2.4593 | 1945.649 | 3092.763 | 3504 |
| channel | global-pruned | drop=3e-02 | 2.0604 | 1944.052 | 3992.833 | 6391 |
| channel | global-pruned | drop=5e-02 | 1.5799 | 1940.802 | 7938.412 | 18888 |
| channel | global-pruned | drop=1e-01 | 0.9935 | 1942.111 | 16325.998 | failed@40000 |
| checker | global-pruned | drop=0e+00 | 96.9969 | 1795.121 | 2287.725 | 228 |
| checker | global-pruned | drop=1e-04 | 7.7556 | 1549.971 | 1647.217 | 231 |
| checker | global-pruned | drop=1e-03 | 5.4953 | 1547.347 | 1643.060 | 234 |
| checker | global-pruned | drop=3e-03 | 4.4253 | 1546.284 | 1630.915 | 235 |
| checker | global-pruned | drop=1e-02 | 3.1966 | 1544.463 | 1810.148 | 815 |
| checker | global-pruned | drop=2e-02 | 2.4822 | 1545.050 | 2550.116 | 3044 |
| checker | global-pruned | drop=3e-02 | 2.0837 | 1544.911 | 3707.782 | 6648 |
| checker | global-pruned | drop=5e-02 | 1.6022 | 1543.214 | 7577.987 | 19649 |
| checker | global-pruned | drop=1e-01 | 1.0059 | 1542.387 | 15910.564 | failed@40000 |

### 4.4 experiment4：支撑选择策略

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | base-local | layers=2 | 6.0922 | 42.157 | 77.549 | 104 |
| continuous | geometric-layer | 2 -> 3 | 12.5471 | 103.474 | 145.762 | 91 |
| continuous | indicator-strong | K=64 one-shot | 6.1082 | 52.314 | 94.384 | 104 |
| continuous | adaptive-budget | R=8, B=128, q=16 | 6.2364 | 150.757 | 189.985 | 101 |
| channel | base-local | layers=2 | 6.0922 | 47.088 | 2837.432 | 6840 |
| channel | geometric-layer | 2 -> 3 | 12.5471 | 115.322 | 3346.021 | 6050 |
| channel | indicator-strong | K=64 one-shot | 6.4855 | 88.460 | 2345.895 | 5592 |
| channel | adaptive-budget | R=8, B=128, q=16 | 6.2390 | 204.632 | 1631.328 | 3553 |
| checker | base-local | layers=2 | 6.0922 | 45.365 | 1679.545 | 4144 |
| checker | geometric-layer | 2 -> 3 | 12.5471 | 133.914 | 277.175 | 243 |
| checker | indicator-strong | K=64 one-shot | 6.4838 | 94.556 | 497.236 | 946 |
| checker | adaptive-budget | R=8, B=128, q=16 | 6.2368 | 194.145 | 531.118 | 807 |

### 4.5 experiment5：有限步全局目标

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| continuous | geometric | P_G | 1.3950 | 2.749 | 49.761 | 131 |
| continuous | Jacobi-global | m=4 | 2.8455 | 16.789 | 68.563 | 112 |
| continuous | Jacobi-global | m=8 | 4.7096 | 38.569 | 82.368 | 100 |
| continuous | Jacobi-global | m=16 | 9.1300 | 116.985 | 168.682 | 91 |
| continuous | Jacobi-global | m=32 | 20.4547 | 508.525 | 573.848 | 86 |
| continuous | Jacobi-global | m=64 | 48.1166 | 2005.076 | 2116.975 | 83 |
| continuous | PCG-global | m=4 | 2.8455 | 44.709 | 73.982 | 90 |
| continuous | PCG-global | m=8 | 4.7096 | 93.887 | 123.564 | 86 |
| continuous | PCG-global | m=16 | 9.1300 | 152.939 | 214.620 | 116 |
| continuous | PCG-global | m=32 | 20.4547 | 331.212 | 393.960 | 93 |
| continuous | PCG-global | m=64 | 48.1166 | 599.376 | 719.806 | 90 |
| continuous | global-exact | tol=1e-10 | 95.7222 | 1602.170 | 1827.152 | 90 |
| channel | geometric | P_G | 1.3950 | 1.222 | 15137.524 | failed@40000 |
| channel | Jacobi-global | m=4 | 2.8455 | 14.235 | 21691.285 | failed@40000 |
| channel | Jacobi-global | m=8 | 4.7096 | 300.130 | 17543.662 | failed@40000 |
| channel | Jacobi-global | m=16 | 9.1300 | 97.705 | 21883.742 | 39142 |
| channel | Jacobi-global | m=32 | 20.4547 | 374.848 | 14466.273 | 18220 |
| channel | Jacobi-global | m=64 | 48.1166 | 1673.239 | 8384.070 | 4704 |
| channel | PCG-global | m=4 | 2.8455 | 42.216 | 8136.200 | 19397 |
| channel | PCG-global | m=8 | 4.7096 | 76.794 | 3790.210 | 9058 |
| channel | PCG-global | m=16 | 9.1300 | 163.501 | 2464.865 | 4043 |
| channel | PCG-global | m=32 | 20.4547 | 266.572 | 537.370 | 347 |
| channel | PCG-global | m=64 | 48.1166 | 528.332 | 1004.407 | 366 |
| channel | global-exact | tol=1e-10 | 97.8420 | 2256.744 | 9304.493 | 3227 |
| checker | geometric | P_G | 1.3950 | 1.697 | 14855.917 | failed@40000 |
| checker | Jacobi-global | m=4 | 2.8455 | 16.561 | 17549.612 | failed@40000 |
| checker | Jacobi-global | m=8 | 4.7096 | 38.942 | 18199.683 | failed@40000 |
| checker | Jacobi-global | m=16 | 9.1300 | 102.816 | 22096.312 | 39538 |
| checker | Jacobi-global | m=32 | 20.4547 | 397.255 | 14657.253 | 18555 |
| checker | Jacobi-global | m=64 | 48.1166 | 1707.928 | 7173.420 | 3974 |
| checker | PCG-global | m=4 | 2.8455 | 42.749 | 8243.533 | 19695 |
| checker | PCG-global | m=8 | 4.7096 | 80.468 | 3437.367 | 7360 |
| checker | PCG-global | m=16 | 9.1300 | 147.681 | 2247.120 | 4040 |
| checker | PCG-global | m=32 | 20.4547 | 347.851 | 623.238 | 347 |
| checker | PCG-global | m=64 | 48.1166 | 694.808 | 1058.916 | 226 |
| checker | global-exact | tol=1e-10 | 96.9969 | 1759.960 | 2356.540 | 228 |

### 4.6 experiment6：通道场 PCG 密集扫描

| Field | Method | Parameter | P density % | Setup ms | Total ms | Cycles |
| --- | --- | --- | --- | --- | --- | --- |
| channel | PCG-global | m=16 | 9.1300 | 181.620 | 2873.227 | 4043 |
| channel | PCG-global | m=18 | 10.3088 | 162.560 | 3463.399 | 3291 |
| channel | PCG-global | m=20 | 11.6236 | 184.544 | 1294.552 | 1861 |
| channel | PCG-global | m=22 | 13.0007 | 218.540 | 934.971 | 1148 |
| channel | PCG-global | m=24 | 14.4398 | 220.927 | 838.038 | 926 |
| channel | PCG-global | m=26 | 15.8038 | 256.735 | 737.212 | 644 |
| channel | PCG-global | m=28 | 17.3034 | 238.089 | 583.501 | 441 |
| channel | PCG-global | m=30 | 18.8537 | 273.400 | 533.431 | 336 |
| channel | PCG-global | m=32 | 20.4547 | 264.775 | 527.054 | 347 |
| channel | PCG-global | m=34 | 21.9609 | 305.936 | 600.780 | 367 |
| channel | PCG-global | m=36 | 23.5999 | 301.991 | 524.300 | 268 |
| channel | PCG-global | m=38 | 25.2786 | 346.370 | 559.545 | 231 |
| channel | PCG-global | m=40 | 26.9969 | 342.108 | 568.008 | 234 |
| channel | PCG-global | m=42 | 28.6040 | 352.662 | 591.283 | 238 |
| channel | PCG-global | m=44 | 30.3389 | 377.080 | 632.433 | 244 |
| channel | PCG-global | m=46 | 32.1029 | 398.410 | 663.274 | 252 |
| channel | PCG-global | m=48 | 33.8960 | 419.569 | 706.833 | 256 |
| channel | PCG-global | m=50 | 35.5642 | 423.573 | 733.491 | 263 |
| channel | PCG-global | m=52 | 37.3532 | 461.425 | 835.872 | 276 |
| channel | PCG-global | m=54 | 39.1612 | 458.902 | 817.514 | 288 |
| channel | PCG-global | m=56 | 40.9881 | 481.286 | 809.324 | 276 |
| channel | PCG-global | m=58 | 42.6795 | 520.740 | 927.034 | 324 |
| channel | PCG-global | m=60 | 44.4826 | 533.322 | 980.107 | 360 |
| channel | PCG-global | m=62 | 46.2950 | 538.426 | 1027.878 | 364 |
| channel | PCG-global | m=64 | 48.1166 | 622.702 | 1172.805 | 366 |
| channel | global-exact | tol=1e-10 | 97.8420 | 2442.930 | 10025.241 | 3227 |

密集扫描显示 PCG 步数与两网格循环数具有明显非单调关系。通道场当前最小循环数位于 $m=38$，为 231 次；$m=36$--50 均低于 300 次，构成较稳定的低循环窗口。$m=30$--34 出现局部反弹，而随后 $m=36$--38 又继续改善，说明自适应停止不能依赖一次相邻检查点变差，必须保留历史最优并设置耐心窗口。与精确全局能量基的 3227 次相比，$m=38$ 的循环数减少约 92.8%，同时 $P$ 密度仅 25.28%。

experiment6 对每个候选都用独立的长双精度稠密 Cholesky 求解确定性粗网格测试系统，并与生产稀疏 Cholesky 结果交叉检查；全部候选通过预设精度阈值。验证工作不计入 `Setup ms` 与 `Total ms`，生产粗求解代码保持不变。因此，在当前测试覆盖范围内，没有证据表明有限步 PCG 的优势来自粗矩阵数值求解错误。

---

## 5. 结果综合与当前结论

### 5.1 主要规律

1. **支撑扩展存在场依赖。** 连续场只需较低复杂度，通道场需要沿长程强连接传播信息，棋盘场更适合统一的中等支撑。
2. **局部求解存在精度平台。** 达到足够精度后，继续迭代不会改善循环数。
3. **全局基具有长尾。** 大量条目可删，但直接构造全局基再剪枝的 setup 代价过高。
4. **全局极小不等于两网格最优。** 固定光滑器下，有限步构造可以形成更匹配的粗空间。
5. **有限步存在宽窗口和局部反弹。** 通道场当前 $m=36$--50 均低于 300 次循环，最低点为 $m=38$ 的 231 次；$m=30$--34 的局部反弹表明单点贪心停止不足以可靠识别窗口。
6. **未发现粗矩阵求解异常。** 独立长双精度交叉检查和稀疏因子压力测试支持“光滑器--粗空间失配”的解释，同时避免把有限个测试右端项夸大为形式证明。
7. **单一固定策略不可能普遍最优。** 自适应方法必须能识别简单场并提前退出，也必须能为少数困难列增加预算。

### 5.2 当前最合理的基线

| 场景 | 当前建议基线 | 原因 |
|---|---|---|
| continuous | 几何插值 | 单次求解总时间最低，额外构造收益有限 |
| channel | PCG 36--40 步 | $m=38$ 循环数最低为 231；$m=36$ 的本次总时间最低 |
| checker | 三层局部能量基 | 243 次循环且 setup 明显低于全局基和高步 PCG |

这些建议只针对当前 $H/h$、对比度和常数右端项，不应固化为最终规则。若目标是可复用预条件器，应优先比较随机子空间估计的近最坏模态因子，而不是直接把当前右端项下的 231 次循环视为普适最优。

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

目标分成两个清晰模式：给定单一右端项时，用当前残差经过 GS 过滤后的实际难消分量评价候选；面向多右端项或通用预条件器时，用固定随机子空间估计近最坏模态。两种模式共享同一套 PCG 检查点、耐心窗口和历史最优回退，但不能混用同一个“最佳步数”结论。

### 6.2 第一阶段：整体 PCG 自适应早停

最先实现一个简单、可解释的整体早停器，不立即进入复杂的逐列决策。以几何插值 $P_0$ 为初值，每次增加 $\Delta m=2$ 或 4 步 PCG，并在检查点 $m_k$ 保存候选 $P_{m_k}$。

在通用预条件器模式下，对少量固定随机向量 $z_\ell$ 先施加一次 GS，使

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

单右端项模式则把 $z_\ell$ 换成当前实际残差及一到两个小扰动探针。两种模式都不使用 PCG 方程残差来判断“是否足够精确”，而直接估计当前粗空间和固定 GS 的互补程度。结合已测得的一次循环成本，预测

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

当前结果最清楚地说明三点：第一，逼近全局能量基的精度与固定 GS 两网格性能并非单调对应；第二，把右端项改为与离散算子无关的常数向量后，通道场仍出现清晰的有限步 PCG 窗口，$m=38$ 为 231 次，而全局精确基为 3227 次，这显著削弱了“现象来自制造右端项”的解释；第三，experiment6 的逐候选长双精度稠密 Cholesky 交叉检查没有发现能够解释该现象的粗矩阵求解错误。下一阶段的核心任务是按使用场景选择实际残差或随机子空间探针，自动识别窗口，并在额外 setup 成本超过预期收益之前返回历史最佳候选。

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
