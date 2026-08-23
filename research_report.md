# 研究方案与阶段报告 v2.3

## 高对比多尺度扩散问题中的稀疏插值与两网格机制

**对应代码：** `two_grids_iteration 2.3.0`  
**研究阶段：** 两网格插值机制研究，不在本版本实现多层递归  
**核心问题：** 在固定粗点集合下，如何用低构造成本和低算子复杂度获得对当前光滑器有效的稀疏插值？

---

## 摘要

考虑高对比、强非均匀扩散方程离散得到的对称正定系统。固定一组粗自由度后，全局能量极小插值具有明确的变分意义，并等于该 C/F 划分下的离散理想插值，但它通常是稠密的。局部能量插值、强连接支撑扩展以及松弛诱导插值都可以看作对全局理想插值的不同稀疏近似。

v2.3 不再扩展大量启发式模块，而把研究收束为三组实验：支撑策略比较、插值构造比较、对比度—尺度稳健性。支撑实验引入“预算匹配的全局剪枝参考”：对于任意候选插值，在全局基中按每列相同的非零元预算选择节点，再重新做能量极小化。它用于区分“预算不足”和“支撑位置错误”。插值实验新增固定粗点约束下的 F 点加权 Jacobi 松弛，以及基于归一化强连接图距离选择支撑、随后局部能量极小化的偏代数方案。

数值结果给出三个主要判断。第一，在 $128/16$ 棋盘场上，local-2 与预算匹配全局剪枝参考使用相同的逐列非零元数量，但循环数分别为 1475 和 217，说明关键问题是支撑位置而非支撑数量；残差预算扩展达到 219 次。第二，Jacobi 松弛从几何插值出发，随步数增加稳定降低 F 点方程残差和全局基能量误差，并在多个单右端项实验中比局部或全局能量基使用更少循环。第三，插值的逐列能量误差并不能单独排序两网格收敛率；两网格性能还取决于插值空间与光滑器未消除误差的匹配。

因此，后续两网格研究不再以“尽可能逼近全局能量基”作为唯一目标，而以三个量联合评价：固定 C/F 约束下的残差—能量误差、两网格误差传播因子、插值与粗算子复杂度。完整多重网格推广推迟到这些关系在两层上得到稳定解释之后。

**关键词：** 高对比扩散；理想插值；能量极小化；Jacobi 松弛；代数距离；稀疏支撑；两网格方法

---

## 1. 问题、范围与应用对象

### 1.1 模型问题

研究

$$
-\nabla\cdot(a(x)\nabla u)=f\quad\text{in }\Omega,
\qquad u=0\quad\text{on }\partial\Omega,
$$

其中

$$
0<a_{\min}\le a(x)\le a_{\max},
\qquad \kappa=a_{\max}/a_{\min}\gg1.
$$

细网格离散得到

$$
Au=f,
$$

$A\in\mathbb R^{n\times n}$ 对称正定。当前代码使用结构网格和五点有限差分，面系数取调和平均；研究结论针对离散 SPD 结构，不能未经验证推广到非对称、非正定问题。

### 1.2 实际系数场

本阶段保留三类具有不同机制的系数场：

| 类型 | 数学特征 | 代表场景 |
|---|---|---|
| 连续随机场 | 对数尺度连续、多频振荡 | 随机介质、热传导参数不确定性 |
| 块状高对比场 | 夹杂与粗尺度块结构 | 复合材料、分区介质 |
| 连通通道场 | 高导结构跨越多个粗单元 | 地下水、油藏渗流、裂隙或高导网络 |

三类名称不是算法路由规则。真正影响局部化的是低能量误差的连通拓扑、粗点约束是否充分，以及这种误差能否被光滑器消除。

### 1.3 研究边界

本版本只研究：

1. 标量 SPD 扩散；
2. 固定几何粗点集合；
3. 一次预平滑、精确 Galerkin 粗解和一次后平滑的两网格框架；
4. 插值支撑、权重与稀疏复杂度。

不同时研究粗点自适应选择、谱增广、多层递归、非结构网格或非对称方程。若全局理想插值仍不能给出可接受的两网格收敛，则应判定固定粗点集合或当前光滑器不足，而不是无限扩大支撑。

---

## 2. 离散理想插值：公式与严格性质

### 2.1 C/F 分块

将未知量重排为 F 点和 C 点：

$$
A=
\begin{bmatrix}
A_{FF}&A_{FC}\\
A_{CF}&A_{CC}
\end{bmatrix}.
$$

由于 $A$ SPD，其主子矩阵 $A_{FF}$ 也 SPD。保持 C 点精确注入的任意插值可写成

$$
P(W)=
\begin{bmatrix}
W\\I
\end{bmatrix},
\qquad W\in\mathbb R^{n_F\times n_C}.
$$

代码中的每个粗基在所属 C 点取 1、其他 C 点取 0，正对应这一形式。

### 2.2 定理 1：全局能量极小插值的唯一表达式

定义总基能量

$$
J(W)=\frac12\operatorname{tr}\bigl(P(W)^TAP(W)\bigr).
$$

则唯一极小点为

$$
W_\star=-A_{FF}^{-1}A_{FC}.
$$

因此

$$
P_\star=
\begin{bmatrix}
-A_{FF}^{-1}A_{FC}\\I
\end{bmatrix}
$$

就是固定 C/F 划分下的理想插值。

**证明。** 展开迹：

$$
J(W)=\frac12\operatorname{tr}(W^TA_{FF}W)
+\operatorname{tr}(W^TA_{FC})
+\frac12\operatorname{tr}(A_{CC}).
$$

对 $W$ 求 Fréchet 导数得到

$$
\nabla J(W)=A_{FF}W+A_{FC}.
$$

驻点满足 $A_{FF}W=-A_{FC}$。$A_{FF}$ SPD，所以 $J$ 严格凸，驻点存在、唯一且为全局极小点。证毕。

### 2.3 定理 2：能量差恒等式

令 $E=W-W_\star$，则

$$
J(W)-J(W_\star)
=\frac12\operatorname{tr}(E^TA_{FF}E).
$$

**证明。** 将 $W=W_\star+E$ 代入 $J$。一次项为

$$
\operatorname{tr}\bigl(E^T(A_{FF}W_\star+A_{FC})\bigr)=0,
$$

剩余二次项即所给表达式。证毕。

该恒等式说明：在相同 C 点注入约束下，候选插值相对全局基的列聚合 A-能量误差，恰好等于它比全局极小值多出的能量，而不是经验指标。

### 2.4 定理 3：支撑受限能量基及嵌套单调性

对第 $j$ 列给定 F 点支撑 $S_j$，定义

$$
w_j^{S_j}
=\arg\min_{\operatorname{supp}(w)\subseteq S_j}
\left(
\frac12 w^TA_{FF}w+w^TA_{FC}e_j
\right).
$$

若 $S_j\subseteq T_j$，则

$$
\|w_j^{T_j}-w_j^\star\|_{A_{FF}}
\le
\|w_j^{S_j}-w_j^\star\|_{A_{FF}}.
$$

**证明。** $S_j$ 对应的可行线性子空间包含于 $T_j$ 对应子空间。扩大支撑后极小目标值不增。由定理 2 的逐列版本，目标值超额等于误差能量的一半，结论成立。证毕。

注意：该定理只保证插值列对全局理想列的能量误差单调，不保证两网格循环数单调。后者依赖整个粗空间与光滑器的组合。

---

## 3. 残差、真实能量误差与支撑选择

### 3.1 定理 4：F 点残差—误差恒等式

定义

$$
R_F(W)=A_{FF}W+A_{FC}.
$$

由 $A_{FF}W_\star+A_{FC}=0$，有

$$
R_F(W)=A_{FF}(W-W_\star)=A_{FF}E.
$$

从而

$$
\operatorname{tr}(E^TA_{FF}E)
=\operatorname{tr}\bigl(R_F(W)^TA_{FF}^{-1}R_F(W)\bigr).
$$

**证明。** 由 $E=A_{FF}^{-1}R_F$ 直接代入左侧。证毕。

这给出了真实全局基误差的精确残差表示，但 $A_{FF}^{-1}$ 本身是全局算子，不能直接作为廉价后验量。

### 3.2 Jacobi 能量残差及其条件界

令 $D=\operatorname{diag}(A_{FF})$，定义

$$
\eta_D(W)^2
=\operatorname{tr}\bigl(R_F(W)^TD^{-1}R_F(W)\bigr)
=\sum_{j}\sum_{i\in F}\frac{r_{ij}^2}{A_{ii}}.
$$

若存在常数 $0<c\le C$ 使

$$
cD\preceq A_{FF}\preceq CD,
$$

则由 SPD 矩阵求逆的序反性，

$$
\frac1C D^{-1}\preceq A_{FF}^{-1}
\preceq\frac1cD^{-1},
$$

进而

$$
\frac1C\eta_D(W)^2
\le \|E\|_{A_{FF},F}^2
\le \frac1c\eta_D(W)^2.
$$

这里的常数可能随对比度和拓扑退化，因此该指标是有条件可靠的误差代理，不是无条件的对比度独立估计。

若 $w_j^{S_j}$ 在支撑上被精确极小化，则一阶最优性给出

$$
(R_F)_i=0,\qquad i\in S_j.
$$

因此它的全部 F 残差位于支撑外。代码的残差预算算法只统计支撑外的 $r_i^2/A_{ii}$；若局部 PCG 未完全收敛，支撑内残差属于代数求解误差，应与局部化误差分开。

### 3.3 当前保留的两种支撑扩展

1. **固定强路径预算。** 从局部基残差和强连接图出发，每列最多加入固定数量的图路径节点。它便宜、可解释，但固定 K 不具备严格尺度一致性。
2. **能量残差预算。** 逐轮标记支撑外 $r_i^2/A_{ii}$，在对称对角归一化的强图上扩展，并只重求解发生变化的列。它构造更贵，但在当前 $128/16$ 棋盘和通道测试中更有效。

旧前沿扩展及旧剪枝评分没有稳定提供新的解释或性能，已从 v2.3 删除。

---

## 4. Jacobi 松弛诱导插值

### 4.1 构造

从任意满足 C 点精确注入的 $W^{(0)}$ 出发，在 F 点执行阻尼 Jacobi：

$$
W^{(m+1)}
=W^{(m)}-\omega D^{-1}
\bigl(A_{FF}W^{(m)}+A_{FC}\bigr).
$$

代码以几何双线性插值为 $W^{(0)}$，每步恢复 C 点单位注入，并在每个 F 行最多保留 8 个绝对值最大的条目。若不截断，这就是对所有理想插值列同时进行的 Jacobi 迭代。

### 4.2 定理 5：未截断 Jacobi 插值的收敛

令

$$
B=D^{-1/2}A_{FF}D^{-1/2},
\qquad
0<\omega<\frac{2}{\lambda_{\max}(B)}.
$$

则

$$
E^{(m+1)}=(I-\omega D^{-1}A_{FF})E^{(m)},
$$

并且

$$
\|E^{(m)}\|_{A_{FF},F}
\le q^m\|E^{(0)}\|_{A_{FF},F},
$$

其中

$$
q=\max_{\lambda\in\sigma(B)}|1-\omega\lambda|<1.
$$

**证明。** 用 $A_{FF}W_\star+A_{FC}=0$ 从迭代式中减去 $W_\star$，得到误差递推。矩阵 $D^{-1}A_{FF}$ 与对称 SPD 矩阵 $B$ 相似，其特征值为正。误差传播矩阵在 $A_{FF}$ 内积下自伴，其谱为 $1-\omega\lambda(B)$；给定阻尼范围使全部特征值绝对值小于 1，由谱定理得到结论。证毕。

### 4.3 定理 6：截断误差传播界

将每步稀疏化造成的扰动记为 $\Delta^{(m)}$：

$$
E^{(m+1)}=TE^{(m)}+\Delta^{(m)},
\qquad T=I-\omega D^{-1}A_{FF}.
$$

若 $\|T\|_{A_{FF}}\le q<1$，则

$$
\|E^{(m)}\|_{A_{FF},F}
\le q^m\|E^{(0)}\|_{A_{FF},F}
+\sum_{\ell=0}^{m-1}q^{m-1-\ell}
\|\Delta^{(\ell)}\|_{A_{FF},F}.
$$

**证明。** 递推展开后应用三角不等式和 $\|T^k\|_{A_{FF}}\le q^k$。证毕。

因此，有限行预算的 Jacobi 插值收敛到理想插值附近的一个误差邻域；要得到严格精度保证，必须控制每步截断的 A-能量，而不能只控制单个条目幅值。

本实现与经典 smoothed aggregation 的共同点是用阻尼松弛改善初始插值；区别是当前 C 点集合预先给定，且只更新 F 行，所以未截断极限恰为 $-A_{FF}^{-1}A_{FC}$。经典 smoothed aggregation 及其收敛分析可参见 Vaněk–Mandel–Brezina 和 Mandel–Brezina 的工作。

---

## 5. 偏代数支撑选择与全局剪枝参考

### 5.1 归一化强连接距离

对矩阵图边 $(i,j)$ 定义

$$
s_{ij}=\frac{|a_{ij}|}{\sqrt{a_{ii}a_{jj}}},
\qquad
\ell_{ij}=\frac{1}{\max(s_{ij},\varepsilon)}.
$$

对每个 F 点计算到 C 点的最短路距离，选择最近的 $q$ 个 C 点；反向汇总得到每个粗列的非凸支撑，随后在该支撑上重新做能量极小化。这样，图距离只负责选择稀疏模式，权重仍由变分问题确定。

该方案比直接逆距离赋权稳定，并在 $64/8$ 通道和棋盘场中以约 3.1 万个非零元获得与 local-2 约 4.4 万非零元相近的循环数。但它不是完整 algebraic-distance AMG：当前没有通过松弛测试向量拟合代数距离，也没有自适应选取 C 点。因此本研究只把它作为“矩阵依赖支撑选择”的低成本比较项，不宣称一般鲁棒性。基于测试向量和最小二乘的更完整方法可参考 Brandt 等人的 algebraic distance 与 bootstrap AMG 工作。

### 5.2 预算匹配的全局剪枝参考

对候选插值 $P$，记第 $j$ 列的 F 点非零元数为 $k_j$。在全局参考列 $p_j^\star$ 中按

$$
\chi_{ij}=|p_{ij}^\star|\sqrt{A_{ii}}
$$

选取最大的 $k_j$ 个 F 点，形成支撑 $S_j^{\mathrm{ref}}$，再在该支撑上重新能量极小化，得到 $P_{\mathrm{ref}(k)}$。

$\chi_{ij}$ 是对角能量贡献代理；由于 $A$ 含非对角耦合，这一 top-k 选择并不严格解出组合意义下的最佳 k 项支撑。因此代码和报告中的 “oracle” 仅表示它使用了不可用于实际 setup 的全局信息，不表示数学上的最优稀疏解。

该参考回答两个不同问题：

1. 若候选与 $P_{\mathrm{ref}(k)}$ 都差，预算可能不足；
2. 若 $P_{\mathrm{ref}(k)}$ 很好而候选差，主要问题是支撑位置或权重选择。

---

## 6. 两网格收敛与插值误差的关系

### 6.1 精确粗校正的变分性质

给定满列秩 $P$，Galerkin 粗校正对应

$$
\Pi_P=P(P^TAP)^{-1}P^TA.
$$

$\Pi_P$ 是到 $\operatorname{range}(P)$ 的 A-正交投影，故

$$
\|(I-\Pi_P)v\|_A
=\min_{y\in\mathbb R^{n_C}}\|v-Py\|_A.
$$

**证明。** 对 $\phi(y)=\|v-Py\|_A^2$ 求导得到正规方程
$P^TAPy=P^TAv$，解为 $y=(P^TAP)^{-1}P^TAv$。残差与 $\operatorname{range}(P)$ A-正交，由勾股恒等式知其为唯一最佳逼近。证毕。

### 6.2 光滑器匹配

设 $S$ 为一次平滑的误差传播算子，$S^\ast$ 为其 A-伴随。对称两网格误差传播为

$$
E_{TG}=S^\ast(I-\Pi_P)S.
$$

若对所有 $w\in\operatorname{range}(P)^{\perp_A}$ 有

$$
\|Sw\|_A\le\mu_P\|w\|_A,
\qquad \mu_P<1,
$$

则

$$
\|E_{TG}\|_A\le\mu_P^2.
$$

**证明。** $I-\Pi_P$ 是 A-正交投影，范数不超过 1；其像属于粗空间的 A-正交补。依次应用 $S$、投影后的假设和 A-伴随的同一范数界即可。证毕。

该结果虽简单，却明确说明：粗空间评价必须针对光滑后仍存在的误差。$P$ 的逐列 A-误差较小，并不自动使 $\mu_P$ 更小。

更标准的两网格理论使用光滑器诱导范数的近似性质，例如

$$
\min_y\|v-Py\|_{\bar R^{-1}}^2
\le K\|v\|_A^2,
$$

并由子空间校正/Xu–Zikatanov 型恒等式得到收敛界。这里 $\bar R$ 是对称化平滑器。重要的是，控制量是“粗空间对平滑误差的逼近”，而不是单独的 $\|P-P_\star\|_{A,F}$。

### 6.3 为什么全局能量基不必给出当前实验的最少循环

全局能量基 $P_\star$ 对固定 C/F 约束最小化列能量，并使 F 点残差为零；它优化的是理想插值/弱逼近性质。当前循环使用全系统 Gauss–Seidel，而不是与这一 C/F 分裂完全匹配的纯 F-relaxation。有限循环数还受右端项中误差模态的权重影响。因此 Jacobi-4 或 Jacobi-8 在某个制造右端项上少于全局参考循环数并不矛盾，也不能据此断言它在算子范数意义下优于全局基。

v2.3 同时输出固定初值的 A-范数误差传播因子代理 `Rho`。它比单一右端项循环数更接近算子性质，但仅使用有限次幂迭代，仍应视为估计而非严格谱半径。

---

## 7. 与 LOD 指数局部化的准确关系

标准 LOD 先给定稳定的局部准插值 $I_H$，定义细尺度核

$$
W=\ker I_H,
$$

再取 W 中的能量正交校正。若准插值具有局部稳定性、局部逼近性，并且核空间满足适当的局部 Poincaré/Caccioppoli 估计，则全局校正满足

$$
\|Qv_H\|_{A,\Omega\setminus U_k}
\le C\theta^k\|Qv_H\|_A,
\qquad 0<\theta<1.
$$

证明的核心是：用逐层截断函数构造核空间测试函数，通过 Galerkin 正交性把外层能量控制为前一层能量的固定比例，递推得到指数衰减。若第 k 层补丁节点数 $K\asymp k^d$，则

$$
\theta^k=\exp(-c k)
=\exp(-cK^{1/d}),
$$

而一般不是对节点数 K 的 $\exp(-cK)$。

当前代码的约束是节点型 C 点注入，细尺度空间是“在 C 点为零”的代数子空间。它与 LOD 都具有“核空间的能量正交补”结构，但并不自动满足标准 LOD 的稳定准插值假设。特别是在高对比连通通道下，局部 Poincaré 常数可能依赖对比度。因此本方案只能把 LOD 衰减作为有条件理论模板，不能把当前实现直接宣称为对比度无关 LOD。

相关基础结果见 Målqvist–Peterseim 的 elliptic multiscale localization；对比度独立局部化需要进一步的加权插值或几何条件。

---

## 8. v2.3 数值实验

### 8.1 脚本划分

| 程序 | 唯一研究方向 |
|---|---|
| `support_comparison` | 所有保留支撑策略及预算匹配全局剪枝参考 |
| `interpolation_comparison` | 几何、Jacobi、强度距离和能量插值 |
| `robustness` | 单种子、单制造右端项下的对比度与尺度变化 |
| `unit_core` | 维数、约束、Jacobi 残差、预算匹配和两网格回归 |

旧九个实验、二进制模型、模型生成器、多右端项、跨种子脚本、README、CHANGELOG 和运行脚本均已移除。系数场在测试程序中确定性生成。

### 8.2 支撑策略：$128/16$、对比度 $10^4$

| 场 | 方法 | nnz(P) | 相对全局能量误差 | 循环 | 同预算全局剪枝循环 |
|---|---|---:|---:|---:|---:|
| 连续 | local-2 | 221089 | 0.0379 | 86 | 86 |
| 连续 | residual-budget | 226321 | 0.0264 | 86 | 86 |
| 通道 | local-2 | 221089 | 0.0586 | >4000 | 2427 |
| 通道 | strong-K64 | 235360 | 0.0359 | 3960 | 2425 |
| 通道 | residual-budget | 226416 | 0.0279 | 2582 | 2426 |
| 通道 | global | 2727127 | 0 | 2416 | — |
| 棋盘 | local-2 | 221089 | 0.0460 | 1475 | 217 |
| 棋盘 | local-3 | 455337 | 0.00852 | 217 | 217 |
| 棋盘 | strong-K64 | 235298 | 0.0261 | 358 | 217 |
| 棋盘 | residual-budget | 226336 | 0.0220 | 219 | 217 |

解释：

- 连续场不存在值得支付的支撑扩展收益。
- 棋盘场中 local-2 与同预算全局剪枝参考之间有巨大差距，残差预算几乎关闭这一差距；其价值主要来自选择正确位置，而不是增加大量非零元。
- 通道场中残差预算接近同预算参考，但同预算参考本身仍约 2426 次，说明固定粗点空间不足是主要剩余瓶颈。
- 全局剪枝参考并非总是优于候选，例如几何预算下可能破坏空间覆盖；因此它只诊断高质量全局列在给定预算下可提供什么，不是通用生产算法。

### 8.3 插值构造：$64/8$、对比度 $10^4$

| 场 | 方法 | nnz(P) | F 残差 | 相对全局能量误差 | 循环 |
|---|---|---:|---:|---:|---:|
| 连续 | geometric | 11025 | $3.97\times10^3$ | 1.314 | 118 |
| 连续 | Jacobi-4 | 21189 | $1.58\times10^3$ | 0.906 | 106 |
| 连续 | Jacobi-8 | 27533 | $1.08\times10^3$ | 0.762 | 95 |
| 连续 | strength-distance | 31409 | $3.02\times10^2$ | 0.139 | 135 |
| 连续 | global | 192129 | $9.00\times10^{-7}$ | 0 | 111 |
| 通道 | geometric | 11025 | $1.67\times10^4$ | 1.313 | 212 |
| 通道 | Jacobi-8 | 27376 | $6.18\times10^3$ | 0.930 | 185 |
| 通道 | strength-distance | 31377 | $6.02\times10^2$ | 0.0500 | 231 |
| 通道 | global | 192129 | $2.85\times10^{-6}$ | 0 | 230 |
| 棋盘 | geometric | 11025 | $1.56\times10^4$ | 1.298 | 200 |
| 棋盘 | Jacobi-8 | 27394 | $6.14\times10^3$ | 0.923 | 148 |
| 棋盘 | strength-distance | 31182 | $5.29\times10^2$ | 0.0368 | 231 |
| 棋盘 | global | 191254 | $3.13\times10^{-6}$ | 0 | 231 |

Jacobi 的残差、全局基误差和循环数随步数总体改善，且 setup 远低于全局能量基。强度距离能以较少非零元逼近全局列，但循环数不随列误差同步改善；它提供了“列误差不等于两网格目标”的直接证据。

### 8.4 对比度与尺度

`robustness` 保留 $64/8$ 和 $128/16$，对比度为 $10^2,10^4,10^6$，固定 seed=1 和一个制造右端项。代表性 $128/16$ 结果为：

| 场/对比度 | geometric | Jacobi-4 | local-3 | strong-K64 | global |
|---|---:|---:|---:|---:|---:|
| 连续/$10^2$ | 66 | 61 | 82 | 82 | 82 |
| 连续/$10^4$ | 87 | 69 | 86 | 86 | 86 |
| 连续/$10^6$ | 146 | 83 | 107 | 106 | 107 |
| 通道/$10^2$ | 686 | 527 | 807 | 965 | 571 |
| 通道/$10^4$ | >4000 | >4000 | >4000 | 3960 | 2416 |
| 通道/$10^6$ | 2040 | 1422 | >4000 | >4000 | 2598 |
| 棋盘/$10^2$ | 505 | 404 | 174 | 174 | 174 |
| 棋盘/$10^4$ | >4000 | >4000 | 217 | 358 | 217 |
| 棋盘/$10^6$ | 378 | 255 | 3395 | >4000 | 218 |

结论不能简化为“更高对比度必然更难”：固定确定性场生成方式下，改变对比度也改变光滑误差的相对权重，非单调结果是可能的。单右端项循环数尤其可能偏向某些模态，故应结合 `Rho` 和全局基误差解释。

---

## 9. 当前价值、创新边界与后续计划

### 9.1 研究价值

本研究的价值不在于提出一个全新的 AMG 算法，而在于建立了一个小而完整的机制验证框架：

1. 全局能量插值、支撑受限插值和 Jacobi 插值由同一个 C/F 理想插值方程统一；
2. 残差与真实全局基误差之间有精确恒等式，Jacobi 指标的适用条件清楚；
3. 预算匹配的全局剪枝参考可以分离支撑数量与支撑位置；
4. 数值结果直接显示列能量误差与两网格收敛目标之间的差异。

### 9.2 创新边界

- 能量最小化插值、理想插值、smoothed aggregation、algebraic distance 和 LOD 都已有成熟研究；这些概念本身不是本项目原创。
- 当前 Jacobi 公式是理想插值方程的标准松弛，不应包装为新算法。
- 可能形成独立研究贡献的是：针对高对比结构，在固定两网格设置下，用预算匹配全局剪枝参考定量诊断支撑选择效率，并据此设计同时面向误差传播与复杂度的自适应准则。该贡献目前仍处于数值证据阶段。

### 9.3 后续只保留两项工作

**工作 A：两网格目标的可计算代理。** 研究

$$
\eta_D(P),\qquad
\|P-P_\star\|_{A,F},\qquad
\widehat\rho(E_{TG}),\qquad
\operatorname{nnz}(P),\operatorname{nnz}(P^TAP)
$$

之间的关系。目标是确定何时继续 Jacobi 松弛、何时转为支撑扩展、何时应停止并判定粗点集合不足。

**工作 B：面向收益/复杂度的单一自适应算法。** 不再增加平行启发式，而从两个有效动作中选择：

1. 对当前插值做少量截断 Jacobi；
2. 对少数高残差列做对称强图支撑扩展和局部能量重优化。

每个动作以

$$
\frac{\Delta\widehat\rho\ \text{或}\ \Delta\eta_D}
{\Delta\operatorname{nnz}(P)+
\gamma\Delta\operatorname{nnz}(P^TAP)}
$$

作为停止和选择依据。只有在 $64/8$、$128/16$ 以及至少一个更大尺度上得到稳定结论后，才讨论多层递归。

---

## 10. 代码结构与复现

代码仅保留：

```text
src/core/linear_algebra.hpp
src/pde/diffusion_problem.hpp
src/multigrid/energy_interpolation.hpp
src/multigrid/algebraic_interpolation.hpp
src/multigrid/adaptive_support.hpp
src/multigrid/residual_budget_support.hpp
src/multigrid/reference_pruning.hpp
src/multigrid/two_grid_solver.hpp
src/experiment/{common,metrics,reporting,test_problem}.hpp
test/{support_comparison,interpolation_comparison,robustness,unit_core}.cpp
```

构建与测试：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

运行：

```bash
./build/support_comparison
./build/interpolation_comparison
./build/robustness
```

三个程序分别生成同名 TXT 和 CSV 结果。`support_comparison` 默认使用 $128/16$；`interpolation_comparison` 默认使用 $64/8$；`robustness --quick` 只运行 $64/8$，不加 `--quick` 同时运行 $64/8$ 与 $128/16$。

---

## 参考文献

1. A. Målqvist, D. Peterseim, “Localization of elliptic multiscale problems,” *Mathematics of Computation*, 83 (2014), 2583–2603. [AMS 原文](https://www.ams.org/journals/mcom/2014-83-290/S0025-5718-2014-02868-8/S0025-5718-2014-02868-8.pdf)
2. P. Vaněk, J. Mandel, M. Brezina, “Algebraic multigrid by smoothed aggregation for second and fourth order elliptic problems,” *Computing*, 56 (1996), 179–196. [期刊页面](https://link.springer.com/article/10.1007/BF02238511)
3. J. Mandel, M. Brezina, “Convergence of algebraic multigrid based on smoothed aggregation,” *Numerische Mathematik*, 88 (2001), 559–579. [技术报告](https://ww3.math.ucla.edu/camreport/cam99-29.pdf)
4. W. L. Wan, T. F. Chan, B. Smith, “An energy-minimizing interpolation for robust multigrid methods,” *SIAM Journal on Scientific Computing*, 21 (2000), 1632–1649. [作者版 PDF](https://cs.uwaterloo.ca/~jwlwan/papers/WanChanSmith00.pdf)
5. J. Brannick, S. P. MacLachlan, J. B. Schroder, B. S. Southworth, “The role of energy minimization in algebraic multigrid interpolation,” 2019. [arXiv:1902.05157](https://arxiv.org/abs/1902.05157)
6. J. Brannick, F. Cao, K. Kahl, R. D. Falgout, X. Hu, “Optimal interpolation and compatible relaxation in classical algebraic multigrid,” *SIAM Journal on Scientific Computing*, 40 (2018), A1473–A1493. [DOI 页面](https://epubs.siam.org/doi/10.1137/17M1123456)
7. A. Brandt, J. Brannick, K. Kahl, I. Livshits, “Bootstrap AMG,” *SIAM Journal on Scientific Computing*, 33 (2011), 612–632. [DOI 页面](https://epubs.siam.org/doi/10.1137/090752973)
8. A. Brandt, J. Brannick, K. Kahl, I. Livshits, “An algebraic distances measure of AMG strength of connection,” 2011. [arXiv:1106.5990](https://arxiv.org/abs/1106.5990)

