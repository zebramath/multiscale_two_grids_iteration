# 研究方案与阶段报告 v2.5

## 全局能量极小插值的空间、代数与稀疏三重局部化

**对应代码：** `two_grids_iteration 2.5.0`  
**研究阶段：** 固定 C/F 划分与对称 Gauss--Seidel 两网格基线  
**统一配置：** $h=1/128$，$H=1/16$，$H/h=8$，对比度 $\kappa=10^4$，系数种子 1  
**本版边界：** 只研究两网格；不引入多层递归、粗点自适应或最终的自适应决策器

---

## 摘要

本研究把同一个全局能量极小插值的三种近似来源放在统一框架中：

1. **空间局部化：** 把每个基函数限制在有限支撑或图邻域；
2. **代数局部化：** 只做有限 PCG/Jacobi 步，或采用有限容差；
3. **稀疏局部化：** 对已构造的全局基进行幅值剪枝。

记 $P_\star$ 为保持粗点精确注入的全局能量极小插值，$S$ 为支撑集合，$m$ 为迭代预算，$\delta$ 为剪枝阈值，则统一记号为

$$
P_{S,m,\delta}.
$$

v2.5 将固定配置统一到 $128/16$，以扩大问题规模并使成本差异更明显；保留原有的 `residual-budget` 支撑扩展；把 Jacobi/PCG 的迭代步数实验分成无剪枝和统一列相对剪枝两组；所有时间改为一次构造、一次两网格运行，不再重复三次取中位数。完整机器可读结果仍保存在代码包的 `results/*.csv`，本文只呈现整理后的汇总表。

本次结果支持以下阶段性判断：

- 支撑从 2 层扩大到 4 层时，基误差可下降一个至两个数量级，但 GS 两网格循环数不一定下降；
- 固定支撑上，容差收紧主要增加 setup，达到误差饱和后几乎不改变循环数；
- 全局基存在显著可剪枝尾部，但先构造全局基的 oracle 成本仍然很高；
- residual-budget 能在相近稀疏度下明显降低 F 残差，且保留了逐列选择和局部重求解的研究接口；
- Jacobi 与 PCG 的“更接近全局基”和“更快两网格”不是同一个排序，剪枝还会改变这一关系。

这些结论是固定网格、固定对比度和单个右端项下的探索性证据，不是对所有系数场或所有多重网格配置的普遍定理。

---

## 1. 模型、两网格和计时口径

考虑

$$
-\nabla\cdot(a(x)\nabla u)=f\quad\text{in }\Omega,
\qquad u=0\quad\text{on }\partial\Omega,
$$

其中 $0<a_{\min}\le a(x)\le a_{\max}$。五点有限差分和面系数调和平均给出 SPD 矩阵 $A$。程序固定三类二维系数场：连续随机场 `continuous`、连通通道场 `channel` 和随机二值棋盘场 `checker`。

所有候选 $P$ 都放入同一个两网格循环：一次前向 GS、一次 Galerkin 粗校正、一次后向 GS。令

$$
\Pi_A(P)=P(P^TAP)^{-1}P^TA,
$$

则误差传播可以写成

$$
E_{\rm TG}(P)=S_b\bigl(I-\Pi_A(P)\bigr)S_f.
$$

因此改变的只有插值，光滑器、粗点集合和停止准则保持不变。

### 1.1 时间定义

$$
T_{\rm setup}=T_P+T_{A_c,\rm factor},
\qquad
T_{\rm total}=T_P+T_{A_c,\rm factor}+T_{\rm solve}.
$$

`Build ms` 是插值构造；`TG setup ms` 包括 $P^TAP$、排序和粗矩阵分解；`Setup ms` 是二者之和；`Total ms` 再加一次右端项求解。每个候选只计时一次。80 步幂迭代只用于估计 $ρ(E_{\rm TG})$，不计入 `Total ms`。循环数来自制造解右端项，属于右端项相关指标；$ρ$ 用于补充渐近行为。

---

## 2. 全局能量极小插值

将未知量重排为 F 点和 C 点：

$$
A=\begin{bmatrix}A_{FF}&A_{FC}\\A_{CF}&A_{CC}\end{bmatrix},
\qquad
P(W)=\begin{bmatrix}W\\I\end{bmatrix}.
$$

由于 $A$ SPD，$A_{FF}$ 也是 SPD。定义基矩阵的总能量

$$
J(W)=\frac12\operatorname{tr}\bigl(P(W)^TAP(W)\bigr).
$$

### 定理 1：全局极小点唯一

$$
W_\star=-A_{FF}^{-1}A_{FC}
$$

是 $J$ 的唯一极小点。

**证明。** 展开迹：

$$
J(W)=\frac12\operatorname{tr}(W^TA_{FF}W)
 +\operatorname{tr}(W^TA_{FC})
 +\frac12\operatorname{tr}(A_{CC}).
$$

其 Fréchet 导数为 $A_{FF}W+A_{FC}$。由于 $A_{FF}$ SPD，$J$ 严格凸；唯一驻点存在且满足 $A_{FF}W+A_{FC}=0$。证毕。

### 定理 2：能量差恒等式

对任意 $W=W_\star+E$，有

$$
J(W)-J(W_\star)=\frac12\operatorname{tr}(E^TA_{FF}E).
$$

**证明。** 将 $W_\star+E$ 代入 $J$。一次项为
$\operatorname{tr}\bigl(E^T(A_{FF}W_\star+A_{FC})\bigr)=0$，剩下二次项即上式。证毕。

代码中的 `global-energy` 用对角预条件 CG 逐列求解到 $10^{-10}$ 相对残差；它是数值参考，不是显式形成逆矩阵的数学精确解。

---

## 3. 三重局部化的理论结构

### 3.1 空间局部化：支撑受限投影

对第 $j$ 个粗基，给定 F 点支撑 $S_j$，令

$$
V_{S_j}=\{v:\operatorname{supp}(v)\subseteq S_j\}.
$$

支撑受限能量基 $w_{S_j}$ 定义为

$$
w_{S_j}=\arg\min_{v\in V_{S_j}}
\left(\frac12v^TA_{FF}v+v^TA_{FC}e_j\right).
$$

其一阶条件为

$$
v^TA_{FF}(w_\star-w_{S_j})=0,
\qquad \forall v\in V_{S_j}.
$$

也就是说，$w_{S_j}$ 是 $w_\star$ 在 $A_{FF}$ 内积下到 $V_{S_j}$ 的正交投影。

### 定理 3：嵌套支撑单调性

若 $S_j\subseteq T_j$，则

$$
\|w_\star-w_{T_j}\|_{A_{FF}}
\le \|w_\star-w_{S_j}\|_{A_{FF}}.
$$

并且对任意 $v\in V_{S_j}$，有

$$
\|w_\star-v\|_{A_{FF}}^2
=\|w_\star-w_{S_j}\|_{A_{FF}}^2
 +\|w_{S_j}-v\|_{A_{FF}}^2.
$$

**证明。** 正交条件给出 $w_\star-w_{S_j}\perp_{A_{FF}}V_{S_j}$。令 $v=w_{S_j}+z$，其中 $z\in V_{S_j}$，展开平方范数即得 Pythagoras 恒等式。取 $v=w_{T_j}$ 并利用 $V_{S_j}\subseteq V_{T_j}$，得到单调性。证毕。

### 3.2 与 LOD 的关系

在 LOD 的局部稳定准插值、核空间稳定分解和局部 Poincaré 条件成立时，corrector 通常满足

$$
\|w_\star-w_k\|_A\le C_{\rm loc}\theta^k,
\qquad 0<\theta<1,
$$

其中 $k$ 是 patch 层数。当前程序的节点注入空间不自动等同于标准 LOD 的 $\ker I_H$，所以这里只把指数衰减作为待验证的条件结论。高对比通道可能使 $C_{\rm loc}$ 和 $\theta$ 依赖通道拓扑；要得到对比度独立结论，需要加权准插值和额外几何假设。

在二维规则网格中，支撑节点数通常满足 $K\asymp k^2$，因此层数指数衰减更自然地写成 $\exp(-c\sqrt K)$，不能未经条件地写成 $\exp(-cK)$。

### 3.3 代数局部化：容差、Jacobi 与 PCG

在固定支撑 $S$ 上，精确局部系统写为

$$
A_{SS}x_S=b_S.
$$

若 $x_{S,m}$ 为有限迭代结果，$r_m=b_S-A_{SS}x_{S,m}$，则有精确恒等式

$$
\|x_S-x_{S,m}\|_{A_{SS}}^2
=r_m^TA_{SS}^{-1}r_m.
$$

因此“容差”和“有限步”只是同一代数误差的两种控制参数。

令 $D=\operatorname{diag}(A_{SS})$，阻尼 Jacobi 为

$$
x^{m+1}=x^m+\omega D^{-1}(b_S-A_{SS}x^m).
$$

由于 $D^{-1}A_{SS}$ 与对称矩阵 $D^{-1/2}A_{SS}D^{-1/2}$ 相似，若

$$
0<\omega<2/\lambda_{\max}(D^{-1}A_{SS}),
$$

则

$$
\|e_m\|_{A_{SS}}\le q_J^m\|e_0\|_{A_{SS}},
\quad
q_J=\max_{\lambda\in\sigma(D^{-1}A_{SS})}|1-\omega\lambda|<1.
$$

对角预条件 CG 满足经典上界

$$
\|e_m\|_{A_{SS}}
\le 2\left(\frac{\sqrt{\kappa_S}-1}{\sqrt{\kappa_S}+1}\right)^m
\|e_0\|_{A_{SS}},
$$

其中 $\kappa_S=\kappa(D^{-1/2}A_{SS}D^{-1/2})$。这是最坏情形上界，不保证欧氏残差或两网格循环数单调下降。

### 引理 4：图传播支撑界

若初始误差支撑在图集合 $G_0$，且预条件器为对角矩阵，则 Jacobi 和 PCG 的第 $m$ 步修正属于

$$
\mathcal K_m(D^{-1}A_{SS},D^{-1}r_0)
$$

生成的空间，从而

$$
\operatorname{supp}(x_m-x_0)\subseteq\mathcal N_m(G_0).
$$

**证明。** 对角缩放不改变非零图；每次乘以 $A_{SS}$ 至多沿稀疏图增加一条边。因此 $A_{SS}^q r_0$ 的支撑包含在 $q$ 次邻域内。Jacobi 是这些向量的递推线性组合，PCG 的每个迭代向量也属于前 $m$ 个 Krylov 向量的线性包，结论由归纳得到。若使用非局部预条件器，该支撑结论失效。证毕。

### 3.4 稀疏局部化：全局基剪枝

对全局列 $p_j$ 定义列相对剪枝

$$
(\mathcal T_\delta p_j)_i=
\begin{cases}
(p_j)_i,& |(p_j)_i|>\delta\max_k|(p_j)_k|,\\
0,&\text{otherwise},
\end{cases}
$$

并始终保留 C 点注入。设 $D_\delta=\mathcal T_\delta(P)-P$，则

$$
\|P_\star-\mathcal T_\delta(P)\|_{A,F}
\le \|P_\star-P\|_{A,F}+\|D_\delta\|_{A,F}.
$$

幅值阈值本身不能给出对比度无关的能量误差界：大量“小”条目可能通过强耦合共同贡献能量。因此代码同时报告全局能量误差、$(AP)_F$ 的对角缩放残差和两网格指标。全局剪枝先构造 $P_\star$，只能作为 oracle 诊断，不能代表低成本 setup。

### 定理 5：三类误差的正交/三角分解

令 $P_S$ 为支撑 $S$ 上的精确极小基，$P_{S,m}$ 为有限迭代基，$\widetilde P_{S,m,\delta}$ 为剪枝结果。则

$$
\|P_\star-P_{S,m}\|_{A,F}^2
=\underbrace{\|P_\star-P_S\|_{A,F}^2}_{E_{\rm space}^2}
 +\underbrace{\|P_S-P_{S,m}\|_{A,F}^2}_{E_{\rm algebra}^2},
$$

以及

$$
\|P_\star-\widetilde P_{S,m,\delta}\|_{A,F}
\le E_{\rm space}+E_{\rm algebra}+E_{\rm drop}.
$$

**证明。** 对每一列，$P_\star-P_S$ 对 $V_S$ 正交，而 $P_S-P_{S,m}\in V_S$，故前两项平方范数按 Pythagoras 相加。剪枝项不再满足同一受限极小化条件，只能对三项使用三角不等式。证毕。

该分解控制的是基误差，并不蕴含两网格谱半径排序。因为 $E_{\rm TG}(P)$ 同时包含光滑器和粗空间投影，一般不存在

$$
\|P_1-P_\star\|_A<\|P_2-P_\star\|_A
\Longrightarrow
\rho(E_{\rm TG}(P_1))<\rho(E_{\rm TG}(P_2)).
$$

这正是有限步 PCG/Jacobi 可能优于更精确全局基的理论原因：有限步改变了粗空间对 GS 难平滑模态的匹配，而非简单地“更接近极小能量”就一定更快。

### 3.5 residual-budget 支撑扩展

`residual_budget_support.hpp` 采用以下可解释流程：

1. 以固定 patch 为初始支撑，计算支撑外的缩放 F 残差；
2. 按列残差平方排序，标记覆盖当前残差能量比例的列；
3. 在强连接图上用优先队列为被标记列添加候选 F 节点；
4. 受每轮节点数、每列总预算和最大轮数限制，只对标记列重新做局部能量最小化；
5. 重算残差，直到达到比例阈值或预算耗尽。

它把“哪里需要扩展”和“扩展多少”分开，并把选择、重求解和残差评估时间都计入 setup。当前实现是研究基线，不声称全局最优；若多数列同时被标记，统一扩层可能更经济。

---

## 4. v2.5 实验矩阵和代码结构

| 脚本 | 只改变的因素 | 主要输出 |
|---|---|---|
| `support_radius.cpp` | 局部层数 $2,3,4$ 与全局 | 稀疏度、基误差、setup/total、循环 |
| `solver_tolerance.cpp` | 固定 3 层上的容差 $10^{-2}$--$10^{-10}$ | 代数误差饱和与成本 |
| `global_pruning.cpp` | 全局列阈值 $0$--$2\times10^{-1}$ | 剪枝稀疏度和两网格退化 |
| `support_strategies.cpp` | local-2、local-3、residual-strong、residual-budget、strength-distance | 支撑选择成本与效果 |
| `iterative_construction.cpp` | Jacobi/PCG 步数 $1,2,4,8,16$，无剪枝与 `drop=1e-2` | 迭代--稀疏--收敛耦合 |

公共实现做了以下效率处理：全局 F 系统只装配一次并复用编号和对角预条件；局部结果使用并行列构造；全局剪枝使用两遍 CSR 扫描；残差预算只重求解被标记列；所有测试共享统一的计时和报告函数。代码包只保留源文件、测试和结果 CSV/TXT，没有冗余 README、CHANGELOG 或旧模型脚本。

---

## 5. 数值结果（$128/16$、$\kappa=10^4$）

表中 `E` 为相对全局能量误差，`P nnz` 和 `dens.` 描述插值稀疏度，`setup/total` 单位为 ms；`failed@4000` 表示达到最大循环数仍未达到 $10^{-6}$。

### 5.1 支撑层数

| 场 | 支撑 | P nnz | dens.% | E | setup | total | cycles |
|---|---:|---:|---:|---:|---:|---:|---:|
| continuous | 2 | 221089 | 6.092 | 3.786e-2 | 57.532 | 98.412 | 86 |
| continuous | 3 | 455337 | 12.547 | 5.231e-3 | 145.954 | 199.239 | 86 |
| continuous | 4 | 745641 | 20.547 | 7.333e-4 | 268.165 | 329.887 | 86 |
| continuous | inf | 3473783 | 95.722 | 0 | 1552.073 | 1732.245 | 86 |
| channel | 2 | 221089 | 6.092 | 5.856e-2 | 61.041 | 2189.975 | failed@4000 |
| channel | 3 | 455337 | 12.547 | 2.018e-2 | 167.023 | 2694.225 | failed@4000 |
| channel | 4 | 745641 | 20.547 | 1.226e-2 | 320.147 | 3165.935 | 3742 |
| channel | inf | 3550710 | 97.842 | 0 | 2144.155 | 7282.485 | 2416 |
| checker | 2 | 221089 | 6.092 | 4.597e-2 | 72.223 | 798.290 | 1475 |
| checker | 3 | 455337 | 12.547 | 8.518e-3 | 186.744 | 323.328 | 217 |
| checker | 4 | 745641 | 20.547 | 1.626e-3 | 332.843 | 504.174 | 217 |
| checker | inf | 3520043 | 96.997 | 0 | 1810.139 | 2367.591 | 217 |

支撑扩大确实降低 $E$，但循环数对场型高度敏感：连续场从 2 层到全局始终为 86 次，棋盘场在 3 层已经降到 217 次，而通道场只有扩大到 4 层或全局后才在 4000 次内收敛。全局基的 96--98% 密度和高 setup 说明“全局精确”不是默认的工程选择。

### 5.2 固定 3 层的 PCG 容差

| 场 | tol | P nnz | E | setup | total | cycles |
|---|---:|---:|---:|---:|---:|---:|
| continuous | 1e-2 | 358176 | 3.352e-2 | 41.477 | 79.336 | 86 |
| continuous | 3e-3 | 429815 | 1.084e-2 | 46.742 | 88.843 | 86 |
| continuous | 1e-3 | 452247 | 6.145e-3 | 53.327 | 90.457 | 86 |
| continuous | 3e-4 | 455303 | 5.326e-3 | 67.248 | 111.234 | 86 |
| continuous | 1e-4 | 455337 | 5.239e-3 | 130.844 | 171.674 | 86 |
| continuous | 1e-10 | 455337 | 5.231e-3 | 135.142 | 174.572 | 86 |
| channel | 1e-2 | 408664 | 3.277e-2 | 83.040 | 2734.376 | failed@4000 |
| channel | 3e-3 | 442840 | 2.162e-2 | 78.806 | 2791.819 | failed@4000 |
| channel | 1e-3 | 450004 | 2.033e-2 | 85.601 | 3102.365 | failed@4000 |
| channel | 3e-4 | 452339 | 2.020e-2 | 100.364 | 2826.802 | failed@4000 |
| channel | 1e-4 | 453779 | 2.019e-2 | 114.512 | 2681.266 | failed@4000 |
| channel | 1e-10 | 455337 | 2.018e-2 | 201.304 | 2757.384 | failed@4000 |
| checker | 1e-2 | 402919 | 2.679e-2 | 51.592 | 395.890 | 583 |
| checker | 3e-3 | 439977 | 1.133e-2 | 62.461 | 192.498 | 221 |
| checker | 1e-3 | 448174 | 8.833e-3 | 71.659 | 228.384 | 217 |
| checker | 3e-4 | 451287 | 8.546e-3 | 94.612 | 210.602 | 217 |
| checker | 1e-4 | 453313 | 8.521e-3 | 84.633 | 217.158 | 217 |
| checker | 1e-10 | 455337 | 8.518e-3 | 180.167 | 310.017 | 217 |

容差收紧后，$P$ 的非零元和能量误差逐步饱和，但循环数几乎不变；在该配置下，$10^{-3}$--$10^{-4}$ 已接近局部化误差主导区。通道场的失败不是容差单独造成的，而是固定 GS 与局部粗空间的组合仍未充分消除通道模态。

### 5.3 全局基相对剪枝（代表阈值）

| 场 | $\delta$ | P nnz | dens.% | E | setup | total | cycles |
|---|---:|---:|---:|---:|---:|---:|---:|
| continuous | 0 | 3473783 | 95.722 | 0 | 1651.193 | 1883.959 | 86 |
| continuous | 1e-3 | 328606 | 9.055 | 1.009e-2 | 1405.372 | 1440.975 | 86 |
| continuous | 1e-2 | 141278 | 3.893 | 7.464e-2 | 1402.661 | 1431.461 | 86 |
| continuous | 5e-2 | 59851 | 1.649 | 3.002e-1 | 1401.412 | 1450.387 | 155 |
| continuous | 2e-1 | 18327 | 0.505 | 9.263e-1 | 1399.597 | 2138.485 | 2030 |
| channel | 0 | 3550710 | 97.842 | 0 | 2141.188 | 7086.871 | 2416 |
| channel | 1e-3 | 204714 | 5.641 | 5.480e-3 | 1890.760 | 2900.938 | 2584 |
| channel | 1e-2 | 115459 | 3.182 | 5.705e-2 | 1889.448 | 2767.990 | 2589 |
| channel | 5e-2 | 57336 | 1.580 | 2.839e-1 | 1887.216 | 3104.602 | failed@4000 |
| channel | 2e-1 | 19357 | 0.533 | 1.002 | 1885.229 | 3407.632 | failed@4000 |
| checker | 0 | 3520043 | 96.997 | 0 | 1895.603 | 2408.135 | 217 |
| checker | 1e-3 | 199426 | 5.495 | 4.505e-3 | 1638.456 | 1723.954 | 218 |
| checker | 1e-2 | 116006 | 3.197 | 5.385e-2 | 1637.608 | 1787.852 | 426 |
| checker | 5e-2 | 58146 | 1.602 | 2.783e-1 | 1647.545 | 2863.267 | failed@4000 |
| checker | 2e-1 | 19657 | 0.542 | 9.977e-1 | 1634.220 | 3165.529 | failed@4000 |

剪枝在 $10^{-3}$--$10^{-2}$ 仍保留可用平台，但通道和棋盘场对较大阈值更敏感。setup 几乎没有随剪枝下降，因为它包含先构造全局基的 oracle 成本；这正是后续“直接构造稀疏局部基”必须解决的问题。

### 5.4 支撑扩展策略

| 场 | 策略 | P nnz | dens.% | E | setup | total | cycles |
|---|---|---:|---:|---:|---:|---:|---:|
| continuous | base-local, 2 | 221089 | 6.092 | 3.786e-2 | 39.002 | 74.020 | 86 |
| continuous | geometric-layer, 2→3 | 455337 | 12.547 | 5.231e-3 | 90.604 | 136.125 | 86 |
| continuous | residual-strong, K=64 | 221669 | 6.108 | 3.780e-2 | 50.581 | 88.278 | 86 |
| continuous | residual-budget, B=128 | 226321 | 6.236 | 2.637e-2 | 156.011 | 192.854 | 86 |
| continuous | strength-distance, q=8 | 127439 | 3.512 | 1.258e-1 | 410.339 | 477.575 | 168 |
| channel | base-local, 2 | 221089 | 6.092 | 5.856e-2 | 43.020 | 2174.399 | failed@4000 |
| channel | geometric-layer, 2→3 | 455337 | 12.547 | 2.018e-2 | 114.829 | 2763.212 | failed@4000 |
| channel | residual-strong, K=64 | 235360 | 6.486 | 3.590e-2 | 88.537 | 2183.139 | 3960 |
| channel | residual-budget, B=128 | 226416 | 6.239 | 2.785e-2 | 157.225 | 1528.761 | 2582 |
| channel | strength-distance, q=8 | 126522 | 3.486 | 5.660e-2 | 560.014 | 2154.367 | failed@4000 |
| checker | base-local, 2 | 221089 | 6.092 | 4.597e-2 | 41.219 | 742.153 | 1475 |
| checker | geometric-layer, 2→3 | 455337 | 12.547 | 8.518e-3 | 102.494 | 241.038 | 217 |
| checker | residual-strong, K=64 | 235298 | 6.484 | 2.606e-2 | 82.966 | 274.887 | 358 |
| checker | residual-budget, B=128 | 226336 | 6.237 | 2.195e-2 | 182.287 | 300.726 | 219 |
| checker | strength-distance, q=8 | 127058 | 3.501 | 5.492e-2 | 555.080 | 2526.418 | failed@4000 |

residual-budget 比 base-local 只增加约 2--3% 的非零元，却在三类场上降低 F 残差；对通道场还把循环数从失败状态降到 2582。代价是额外的残差、图搜索和局部重求解，连续场和棋盘场的总时间仍不如直接 local-2。它因此适合作为“仅对异常列投入预算”的后续自适应原型，而不是无条件替代固定层数。

### 5.5 Jacobi/PCG 固定步数：无剪枝

下表选取 $m=4,8,16$ 展示趋势；$m=0$ 的几何初值以及完整 $m=1,2$ 数据见 `iterative_construction.csv`。

| 场 | 方法 | m | P nnz | E | setup | total | cycles |
|---|---|---:|---:|---:|---:|---:|---:|
| continuous | Jacobi | 4 | 103265 | 8.348e-1 | 17.969 | 48.179 | 69 |
| continuous | Jacobi | 8 | 170913 | 7.002e-1 | 43.856 | 74.290 | 65 |
| continuous | Jacobi | 16 | 331329 | 5.378e-1 | 110.165 | 145.743 | 63 |
| continuous | PCG | 4 | 103265 | 4.836e-1 | 48.829 | 76.812 | 61 |
| continuous | PCG | 8 | 170913 | 2.812e-1 | 79.533 | 112.621 | 66 |
| continuous | PCG | 16 | 331329 | 1.024e-1 | 129.175 | 190.864 | 103 |
| channel | Jacobi | 4 | 103265 | 1.094 | 15.348 | 1710.914 | failed@4000 |
| channel | Jacobi | 8 | 170913 | 9.650e-1 | 38.314 | 1654.889 | failed@4000 |
| channel | Jacobi | 16 | 331329 | 7.896e-1 | 112.456 | 2097.678 | failed@4000 |
| channel | PCG | 4 | 103265 | 6.994e-1 | 41.668 | 1864.442 | failed@4000 |
| channel | PCG | 8 | 170913 | 4.453e-1 | 70.137 | 1793.638 | 3734 |
| channel | PCG | 16 | 331329 | 2.000e-1 | 128.478 | 759.200 | 1125 |
| checker | Jacobi | 4 | 103265 | 1.086 | 9.303 | 1618.608 | failed@4000 |
| checker | Jacobi | 8 | 170913 | 9.570e-1 | 36.399 | 1847.459 | failed@4000 |
| checker | Jacobi | 16 | 331329 | 7.815e-1 | 93.818 | 2113.846 | failed@4000 |
| checker | PCG | 4 | 103265 | 6.858e-1 | 41.465 | 1929.593 | failed@4000 |
| checker | PCG | 8 | 170913 | 4.380e-1 | 89.544 | 1489.500 | 3184 |
| checker | PCG | 16 | 331329 | 1.906e-1 | 136.263 | 705.063 | 962 |

PCG 在相同步数下通常更快降低能量误差，但最佳总时间和最佳循环数不一致：连续场的早期 Jacobi/PCG 已足够，通道和棋盘场需要更深的 PCG 传播才能进入可收敛区。全局参考虽误差为零，却分别需要 1730、7303 和 2125 ms 左右的总时间（具体单次数据在 CSV 中）。

### 5.6 迭代后统一剪枝的影响（$m=8$）

这里对 Jacobi 和 PCG 的同一步无剪枝结果都施加相同的列相对阈值 `drop=1e-2`；剪枝时间计入 setup。

| 场 | 方法 | 版本 | P nnz | E | setup | total | cycles |
|---|---|---|---:|---:|---:|---:|---:|
| continuous | Jacobi | none | 170913 | 7.002e-1 | 43.856 | 74.290 | 65 |
| continuous | Jacobi | drop | 70156 | 7.056e-1 | 43.334 | 72.672 | 65 |
| continuous | PCG | none | 170913 | 2.812e-1 | 79.533 | 112.621 | 66 |
| continuous | PCG | drop | 110945 | 2.936e-1 | 80.453 | 109.505 | 66 |
| channel | Jacobi | none | 170913 | 9.650e-1 | 38.314 | 1654.889 | failed@4000 |
| channel | Jacobi | drop | 68840 | 9.680e-1 | 36.620 | 1469.078 | failed@4000 |
| channel | PCG | none | 170913 | 4.453e-1 | 70.137 | 1793.638 | 3734 |
| channel | PCG | drop | 103989 | 4.534e-1 | 70.926 | 1785.427 | failed@4000 |
| checker | Jacobi | none | 170913 | 9.570e-1 | 36.399 | 1847.459 | failed@4000 |
| checker | Jacobi | drop | 68846 | 9.598e-1 | 35.314 | 1635.385 | failed@4000 |
| checker | PCG | none | 170913 | 4.380e-1 | 89.544 | 1489.500 | 3184 |
| checker | PCG | drop | 104363 | 4.456e-1 | 91.731 | 1559.680 | 3427 |

在连续场中剪枝几乎不改变循环数；通道场的 PCG-$m=8$ 剪枝后略低于 $4000$ 次上限，说明稀疏扰动可能把边界状态推入失败区。棋盘场则显示 P nnz 减少约 39% 后循环数增加。剪枝影响必须和系数拓扑、迭代深度一起测量，不能只看非零元数量。

---

## 6. 结果解释、研究价值和创新边界

### 6.1 当前研究重点

当前最清晰、最小而完整的主线是：

> 在固定两网格和固定 GS 光滑器下，用空间支撑、代数迭代和稀疏剪枝三种局部化共同逼近全局能量极小插值，并以“基误差—稀疏度—setup—两网格循环”四元指标识别主导误差来源。

这条主线同时解释了为什么全局精确基通常不是工程最优，也解释了为什么局部 PCG/Jacobi 的中间步数可能优于严格全局求解。

### 6.2 价值

- **方法学价值：** 将支撑、容差/迭代和剪枝拆成可重复的单变量实验，而不是把误差变化混在一个算法里；
- **工程价值：** 直接报告 setup 和单右端项总时间，能判断额外局部化是否值得；
- **理论价值：** 支撑投影正交性、Krylov 图传播和剪枝三角界把 LOD 风格的空间指数衰减与 Jacobi/PCG 的迭代衰减放进同一误差框架；
- **后续价值：** residual-budget 提供了按列分配支撑扩展和局部求解预算的最小原型。

### 6.3 创新性边界

能量最小插值、Krylov 构造、LOD 局部化和幅值剪枝本身都已有成熟研究。单纯把它们组合在一起，创新性属于有限的框架性探索。更有潜力的后续创新点是：

1. 用可计算的支撑外/支撑内残差区分 $E_{\rm space}$ 与 $E_{\rm algebra}$；
2. 用 residual-budget 或其廉价代理在 setup 预算内逐列分配扩展、迭代和剪枝；
3. 解释固定 GS 下“最佳有限迭代窗口”的可观测代理，而不是只追求 $P_\star$；
4. 在 $h/H$、对比度和通道拓扑变化下验证误差和成本是否具有稳健标度；
5. 将两网格结论递归到多层，并检查总复杂度和误差预算是否保持。

---

## 7. 后续研究计划（保持两网格主线）

1. **先做稳健性：** 只增加 $h/H$、对比度和少量种子，不增加新的插值家族；
2. **再做廉价误差识别：** 每列只计算一次缩放 F 残差和候选边界残差，设定“扩展/继续迭代/剪枝/停止”门槛；
3. **再做成本模型：** 用已测 `Build ms`、粗矩阵 nnz 和每循环应用成本估计剩余总时间；
4. **最后再推广多层：** 先在每个层级复用同一两网格决策器，再分析层间误差和复杂度，不在当前版本同时引入新的光滑器或粗点策略。

一个可检验的自适应目标是

$$
\min_{S,m,\delta}
T_{\rm setup}(S,m,\delta)+T_{\rm solve}(S,m,\delta)
$$

并满足基误差或两网格循环的约束；自适应选择本身的残差、图搜索和重求解时间必须计入 $T_{\rm setup}$。当被标记列比例过高时，应退化为统一扩层；当剪枝扰动超过误差预算时，应停止剪枝；当三类误差都很小而循环仍差时，应检查 C 点集合或光滑器，而不是继续扩大同一插值。

---

## 8. 参考文献

1. W. L. Wan, T. F. Chan, B. Smith, “An Energy-Minimizing Interpolation for Robust Multigrid Methods,” *SIAM Journal on Scientific Computing*, 2000. <https://doi.org/10.1137/S1064827598334277>
2. P. S. Vassilevski, “General Constrained Energy Minimization Interpolation Mappings for AMG,” *SIAM Journal on Scientific Computing*, 2010. <https://doi.org/10.1137/080726252>
3. L. N. Olson, J. B. Schroder, R. S. Tuminaro, “A General Interpolation Strategy for Algebraic Multigrid Using Energy Minimization,” *SIAM Journal on Scientific Computing*, 2011. <https://doi.org/10.1137/100803031>
4. A. Målqvist, D. Peterseim, “Localization of Elliptic Multiscale Problems,” *Mathematics of Computation*, 2014. <https://doi.org/10.1090/S0025-5718-2014-02868-8>
5. F. Hellman, A. Målqvist, “Contrast Independent Localization of Multiscale Problems,” *Multiscale Modeling & Simulation*, 2018. <https://doi.org/10.1137/16M1100460>
6. C. Janna, A. Franceschini, J. B. Schroder, L. Olson, “Parallel Energy-Minimization Prolongation for Algebraic Multigrid,” *SIAM Journal on Scientific Computing*, 2023. <https://doi.org/10.1137/22M1513794>

---

## 附：复现说明

在代码目录中直接编译所有 `test/*.cpp`（C++17，`-O3`，可选 OpenMP），然后运行：

```bash
./build-direct/unit_core
./build-direct/support_radius --threads=4
./build-direct/solver_tolerance --threads=4
./build-direct/global_pruning --threads=4
./build-direct/support_strategies --threads=4
./build-direct/iterative_construction --threads=4
```

默认配置已经是 `--fine=128 --coarse=16 --contrast=1e4`；命令行参数仍可用于小规模调试。完整行级数据位于 `results/support_radius.csv`、`solver_tolerance.csv`、`global_pruning.csv`、`support_strategies.csv` 和 `iterative_construction.csv`。
