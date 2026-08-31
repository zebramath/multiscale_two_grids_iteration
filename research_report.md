# 多尺度有限 Krylov 能量插值研究方案 v6.1

## 1. 研究定位

本项目研究高对比扩散问题中，全局能量极小插值的有限 PCG 路径与实际两网格收敛
目标之间的非单调错位。中心现象是：插值能量随 PCG 迭代持续下降，而固定磨光器下的
两网格循环数可以先快速下降、在有限步达到最小、随后回升，并在高精度能量终点显著
增大。

研究集中回答三个问题：

1. 能量距离、Galerkin 粗矩阵、粗空间几何和两网格谱因子之间有什么严格联系；
2. 有利有限 Krylov 误差方向何时存在，以及性能收益与能量代价分别是什么阶；
3. 最优有限步在什么条件下具有 $\Theta(h^{-1})$ 尺度，以及如何据此构造低 setup 的
   尺度自适应规则。

二维规则网格、固定粗点、对称 Gauss--Seidel、精确 Galerkin 粗校正和两网格结构构成
受控研究平台，使有限 Krylov 插值对粗空间质量的作用能够被单独识别。

### 1.1 相关工作与本文差异

能量极小插值从能量、有约束近核空间再现和稀疏支撑构造延伸到并行实现，形成了成熟的
AMG 插值路线（[Wan--Chan--Smith 2000](https://doi.org/10.1137/S1064827598334277)、
[Olson--Schroder--Tuminaro 2011](https://doi.org/10.1137/100803031)、
[Janna 等 2023](https://doi.org/10.1137/22M1513794)）。另一方面，直接以两网格收敛率
定义最优插值说明能量代理与求解目标应明确区分
（[Brannick 等 2018](https://doi.org/10.1137/17M1123456)）。
JSC 中关于谱粗空间的研究也表明，针对算子范数构造的标准粗空间未必最小化两层迭代的
谱半径，这进一步说明代理目标与实际渐近收敛目标需要区分
（[Ciaramella--Vanzan 2022](https://doi.org/10.1007/s10915-022-01840-9)）。

本项目研究同一能量方程的有限 Krylov 路径：先刻画能量终点与两网格目标的结构性
错位，再由最优步数的尺度估计构造单次有限步插值。

## 2. 模型问题与比较对象

考虑

$$
-\nabla\cdot(a(x)\nabla u)=f\quad\text{in }\Omega=(0,1)^2,
\qquad u=0\quad\text{on }\partial\Omega,
$$

其中 $a(x)>0$，对比度取 $10^2$、$10^4$ 或 $10^6$。节点型有限差分离散得到
对称正定系统 $A_hu_h=b_h$。系数结构包括交叉、曲折、对角、平行、分支和闭合曲环
六类高导通道，并叠加由 seed 控制的块状高导背景。
系数几何以网格单元数定义，因此尺寸轴表示一组受控离散问题，而不是保持同一连续
系数场物理几何不变的标准网格细化实验。

按细点和粗点分块，写成

$$
A=\begin{bmatrix}B&C\\ C^T&A_{CC}\end{bmatrix},
\qquad
P(W)=\begin{bmatrix}W\\I\end{bmatrix}.
$$

比较三种方法：

1. `adaptive`：尺度自适应的有限 PCG 插值；
2. `global-reference`：全局能量方程达到相对欧氏残量 $10^{-10}$ 的高精度数值参考；
3. `geometric`：规则粗点上的双线性几何插值 $P_G$。

理论中的精确能量极小点记为 $W_*$，`global-reference` 是它的高精度数值实现。

## 3. 有限 Krylov 插值机制

能量泛函为

$$
J(W)=\frac12\operatorname{tr}\bigl(P(W)^TAP(W)\bigr).
$$

其唯一极小点与 Schur 补为

$$
W_*=-B^{-1}C,
\qquad
S=A_{CC}-C^TB^{-1}C,
$$

并有精确恒等式

$$
J(W)-J(W_*)=\frac12\lVert W-W_*\rVert_{B,F}^2,
\qquad
P(W)^TAP(W)=S+(W-W_*)^TB(W-W_*).
$$

PCG 路径 $W_m$ 不仅在 $B$ 能量范数下单调逼近 $W_*$，而且任意两个路径位置
$m<\ell$ 之间具有精确尾和

$$
J(W_m)-J(W_\ell)=\frac12\sum_{k=m}^{\ell-1}
\lVert W_{k+1}-W_k\rVert_{B,F}^2.
$$

因此每段插值能量下降都能逐步归因到实际 PCG 更新，而不是只由一个渐近界控制。
两网格误差传播算子

$$
E_{TG}(W)=G^{\dagger_A}
\left[I-P(W)(P(W)^TAP(W))^{-1}P(W)^TA\right]G
$$

同时取决于粗空间相对磨光后误差的方向。能量衡量总体距离，两网格收敛由最坏可见
模态的捕获决定，这一目标差异产生有限步插值优势。

令

$$
Z=B^{1/2}(W-W_*)S^{-1/2}.
$$

则

$$
2\{J(W)-J(W_*)\}=\lVert ZS^{1/2}\rVert_F^2,
$$

$$
P(W)^TAP(W)=S^{1/2}(I+Z^TZ)S^{1/2},
$$

且候选粗空间与能量终点粗空间的主角满足

$$
\tan\theta_i=\sigma_i(Z).
$$

归一化 Galerkin 粗矩阵为 $I+Z^TZ\succeq I$，其非平凡特征值为
$1+\sigma_i(Z)^2$，行列式膨胀因子为 $\prod_i(1+\sigma_i(Z)^2)$。图空间投影及其
正交补还可同时闭式写成

$$
\Pi(Z)=\begin{bmatrix}Z\\I\end{bmatrix}(I+Z^TZ)^{-1}
       \begin{bmatrix}Z^T&I\end{bmatrix},
\qquad
I-\Pi(Z)=N(Z)N(Z)^T,
$$

$$
N(Z)=\begin{bmatrix}I\\-Z^T\end{bmatrix}(I+ZZ^T)^{-1/2}.
$$

在欧氏能量坐标中，若 $T=A^{1/2}GA^{-1/2}$，则

$$
\rho_{TG}(W)=\lVert(I-\Pi(Z))T\rVert_2^2.
$$

因此 $Z$ 同时连接 PCG 能量尾、Galerkin 粗矩阵、粗空间主角和两网格投影目标。终点
附近的一阶结论允许主导奇异值具有任意有限重数：若 $U,V$ 是主导左右奇异子空间，
$U_F$ 和 $Y$ 是其细空间可见性与粗空间响应坐标，则方向 $H$ 上有

$$
\left.\frac{d}{dt}\sqrt{\rho_{TG}(tH)}\right|_{t=0^+}
=-\lambda_{\min}\!\left(\operatorname{sym}(U_F^THY)\right).
$$

因此有限步改善条件作用于整个主导子空间，不依赖非唯一奇异向量基的选择。更强地，
由于主导左奇异模态位于能量终点粗空间的正交补，其细坐标自动满足
$U_F^TU_F=I_r$。严格有利方向存在当且仅当

$$
\operatorname{rank}(Y)=r.
$$

这个唯一的满列秩条件表示主导模态的粗空间响应彼此独立，不再需要额外的细空间
可见性假设。条件成立时可显式取

$$
H_0=U_F(Y^TY)^{-1}Y^T,
\qquad U_F^TH_0Y=I_r,
$$

于是方向导数严格等于 $-1$。这个 $H_0$ 还是所有满足 $U_F^THY=I_r$ 的方向中唯一的
Frobenius 最小范数解，且
$\lVert H_0\rVert_F^2=\operatorname{tr}((Y^TY)^{-1})$。标量公式是重数为一时的特例。

若 $g(Z)=\lambda_{\min}(\operatorname{sym}(U_F^TZY))$，则终点附近存在统一常数
$K,r_0>0$，使

$$
\sqrt{\rho_{TG}(Z)}=\sqrt{\rho_{TG}(0)}-g(Z)+R(Z),
\qquad |R(Z)|\le K\lVert Z\rVert_F^2.
$$

因此在满足 $g(Z)\ge\mu\lVert Z\rVert_F$ 的有利锥中，只要
$\lVert Z\rVert_F\le\min\{r_0,\mu/(2K)\}$，便有

$$
\sqrt{\rho_{TG}(0)}-\sqrt{\rho_{TG}(Z)}
\ge\frac\mu2\lVert Z\rVert_F,
$$

而

$$
\frac{\lambda_{\min}(S)}2\lVert Z\rVert_F^2
\le J(W)-J(W_*)
\le\frac{\lambda_{\max}(S)}2\lVert Z\rVert_F^2.
$$

这把中心现象定量化为“性能改善是一阶量，能量超额是二阶量”。逆向锥
$g(Z)\le-\mu\lVert Z\rVert_F$ 中则得到同阶的性能恶化下界。对逐步收缩且方向相干的
PCG 局部尾，线性项最终严格压过二阶余项，从而严格推出“能量下降、两网格因子上升”。

全局稳定界还能把固定性能优势变成定量能量屏障。若候选相对能量终点的
$\sqrt{\rho_{TG}}$ 改善至少为 $\gamma$，则

$$
J(W)-J(W_*)\ge
\frac{\lambda_{\min}(S)}2
\frac{\gamma^2}{\lVert T\rVert_2^2-\gamma^2}.
$$

因此具有固定幅度优势的有限步候选不能无限逼近能量终点。

反方向也有显式传递。记 $\delta_J=J(W)-J(W_*)$，则

$$
\left|\sqrt{\rho_{TG}(W)}-\sqrt{\rho_{TG}(W_*)}\right|
\le \lVert T\rVert_2
\sqrt{\frac{2\delta_J}{\lambda_{\min}(S)+2\delta_J}}.
$$

只要该上界与终点的 $\sqrt{\rho_{TG}}$ 之和小于 1，就同时给出候选两网格收敛性和
达到指定能量误差比例所需循环数的显式证书。于是全局理论形成双向闭环：能量接近度
控制性能偏差，而固定性能优势反过来要求非零能量距离。

## 4. 尺度自适应选择

令 $n=1/h$ 为每个方向的细网格区间数，并定义矩阵对角尺度比

$$
\chi_A=\frac{\max_i A_{ii}}{\min_i A_{ii}}.
$$

### 4.1 最优步数的条件性尺度估计

令 $D_h=\operatorname{diag}(B_h)$，

$$
\kappa_h=\operatorname{cond}_2(D_h^{-1/2}B_hD_h^{-1/2}),
\qquad
d_{m,h}=\frac{\lVert W_{m,h}-W_h^*\rVert_{B_h,F}}
{\lVert W_{0,h}-W_h^*\rVert_{B_h,F}}.
$$

标准 PCG 上界为

$$
d_{m,h}\le 2\left(
\frac{\sqrt{\kappa_h}-1}{\sqrt{\kappa_h}+1}
\right)^m.
$$

若 $\kappa_h\le\kappa_+h^{-2}$，把能量误差降低到任意固定比例 $\eta$ 只需

$$
m\ge \frac12\left(\frac{\sqrt{\kappa_+}}h+1\right)
\log\frac2\eta=O(h^{-1}).
$$

进一步，只假设相关 PCG 谱分量不发生随细化消失的超快衰减，即在最优点以前

$$
d_{m,h}\ge a_-e^{-\gamma_+mh},
$$

并且两网格最优点处于与 $h$ 无关的可见误差带
$\eta_-\le d_{m_{\mathrm{opt}}(h),h}\le\eta_+$，则

$$
\frac1{\gamma_+}\log\frac{a_-}{\eta_+}\,h^{-1}
\le m_{\mathrm{opt}}(h)\le
\frac12\left(\frac{\sqrt{\kappa_+}}h+1\right)
\log\frac2{\eta_-}.
$$

这给出 $m_{\mathrm{opt}}(h)=\Theta(h^{-1})$ 的条件性估计：标准 PCG 条件数界给出上界，
路径非退化条件给出下界。常数允许随固定的对比度范围、拓扑类别和粗化比变化；PCG 的
代数终止上界仍由细点系统维数决定。七个离线 step--2 sampled oracle 的
$m_{\mathrm{or}}/(1/h)$ 位于
$0.125$--$0.484$，支持当前 $[1/8,1/2]$ 归一化窗口；这些离散采样点是尺度证据，
不等同于连续候选路径上的 $m_{\mathrm{opt}}(h)$。

更具体地，上述双侧界还给出归一化最优步的渐近夹逼

$$
\frac1{\gamma_+}\log\frac{a_-}{\eta_+}
\le\liminf_{h\to0}hm_{\mathrm{opt}}(h)
\le\limsup_{h\to0}hm_{\mathrm{opt}}(h)
\le\frac{\sqrt{\kappa_+}}2\log\frac2{\eta_-}.
$$

### 4.2 尺度自适应规则

令 $n_H=1/H$。adaptive 直接取

$$
m_{\mathrm{ad}}=
\begin{cases}
\operatorname{round}(n/8),&n_H\le8,\\
\operatorname{round}(n/4),&n_H>8,\ \chi_A<10^3,\\
\operatorname{round}(n/3),&n_H>8,\ 10^3\le\chi_A<10^5,\\
\operatorname{round}(n/2),&n_H>8,\ \chi_A\ge10^5.
\end{cases}
$$

该规则只使用粗空间分辨率和矩阵对角尺度比；$10^3$ 与 $10^5$ 是三个测试对比度的
对数中点。每个分支都满足 $m_{\mathrm{ad}}=\operatorname{round}(\theta/h)$，
$\theta\in\{1/8,1/4,1/3,1/2\}$，因此与 $\Theta(h^{-1})$ 尺度窗口一致，整数舍入的
归一化误差不超过 $h/2$。算法只推进一次 PCG 路径、组装一次插值并建立一次粗算子。

## 5. 实验设计与指标

所有方法使用相同矩阵、粗点、右端项、前后 Gauss--Seidel 和 Galerkin 粗算子。正式
求解从零初值开始，相对欧氏残量容差为 $10^{-6}$。adaptive 和 global-reference 的
循环上限为 20000；为更充分地区分缓慢收敛和未收敛，geometric 的上限为 30000。

核心证据包括：

1. 13 个问题比较尺寸、对比度和通道拓扑，其中两组 256/16 作为大尺度扩展；
2. 3 条归一化有限 PCG 路径比较能量超额与两网格循环数；
3. 7 个设计问题和 3 个参数冻结后验证问题上的 step--2 窗口受限离线 oracle；
4. 5 个系数 seed、6 个迁移 RHS 和中心 128/16 问题的 5 次重复计时。

主要指标为循环数、收敛状态、插值密度和以实际循环数 $k$ 定义的有效收敛因子

$$
\rho_{\mathrm{eff}}=r_k^{1/k},
$$

结果文件还给出最后 32 次循环的尾部因子。达到上限且残量仍收缩的结果记为
`slow-limit`；残量
非有限或持续增长的结果记为 `diverged`。时间只在中心 128/16 问题上比较两个完整
收敛的方法；计时入口仅在二者均达到正式容差后生成表格，并报告预热后五次测量的
Q1、中位数和 Q3。

## 6. 跨尺寸、对比度与拓扑的主结果

| 轴 | $1/h,1/H$ | 对比度/拓扑 | adaptive `m/cycles` | reference cycles | geometric |
|---|---|---|---:|---:|---:|
| size | 32,8 | $10^4$/cross | 4 / 37 | 54 | 40 |
| size | 64,8 | $10^4$/cross | 8 / 123 | 264 | 302 |
| size | 64,16 | $10^4$/cross | 21 / 65 | 321 | 24950 |
| center | 128,16 | $10^4$/cross | 43 / 242 | 3227 | 30000 slow |
| large | 256,16 | $10^4$/cross | 85 / 998 | 13948 | 30000 slow |
| large | 256,16 | $10^4$/ring | 85 / 1002 | 3950 | 30000 slow |
| contrast | 128,16 | $10^2$/cross | 32 / 293 | 756 | 1501 |
| contrast | 128,16 | $10^6$/cross | 64 / 273 | 3471 | 30000 slow |
| topology | 128,16 | $10^4$/meandering | 43 / 266 | 759 | 30000 slow |
| topology | 128,16 | $10^4$/diagonal | 43 / 346 | 418 | 30000 slow |
| topology | 128,16 | $10^4$/parallel | 43 / 260 | 421 | 30000 slow |
| topology | 128,16 | $10^4$/branching | 43 / 317 | 615 | 30000 slow |
| topology | 128,16 | $10^4$/ring | 43 / 237 | 768 | 30000 slow |

adaptive 和 global-reference 均在 13 个问题上收敛，累计循环数分别为 4459 和 28972；
global-reference 为 adaptive 的 6.50 倍。geometric 在 4 个问题上收敛，另外 9 个在
30000 次时仍为 `slow-limit`，没有发散结果。

尺度伸缩在 256/16 上保持有效：cross-channel 的 adaptive 为 998 次，而 reference
为 13948 次；winding-ring 为 1002 次和 3950 次。中心问题的有效
收敛因子为 0.944407，对应 reference 的 0.995727；256 cross 上两者为 0.986251 和
0.999010。收敛率指标与循环数给出一致判断。

## 7. 有限路径证据

| 问题 | 采样最优点 | 有限循环 | reference 循环 | 归一化能量超额 |
|---|---:|---:|---:|---:|
| 64/16, cross, $10^4$ | $m=16$ | 60 | 321 | $7.372\times10^{-4}$ |
| 128/16, cross, $10^4$ | $m=40$ | 234 | 3227 | $3.030\times10^{-4}$ |
| 128/16, cross, $10^6$ | $m=64$ | 273 | 3471 | $4.656\times10^{-5}$ |

三条路径均使用 $m/n=0,1/16,\ldots,1/2$ 的归一化检查点。能量超额沿路径持续下降，
循环数则在有限区间达到较小值后回升。128/16、$10^4$ 问题从 $m=40$ 的 234 次上升
到 $m=64$ 的 366 次，最终 reference 为 3227 次；64/16 问题从 $m=16$ 的 60 次
上升到 $m=32$ 的 114 次，最终 reference 为 321 次。这构成目标错位的直接数值证据。

## 8. 窗口受限离线 sampled oracle 选择质量

离线 sampled oracle 检查 $m=0$ 以及从 $n/8$ 到 $n/2$ 的 step--2 路径点。七个设计
问题用于确定归一化区间，三个参数冻结后的新组合用于独立检查；oracle 仅指预设窗口内
该离散候选集合的最优收敛点，不表示完整 PCG 路径上的全局最优点。

| 问题 | adaptive `m/cycles/gap` | oracle `m/(1/h)/cycles` |
|---|---:|---:|
| 32/8, $10^4$, cross | 4 / 37 / 0.00% | 4 / 0.125 / 37 |
| 64/8, $10^4$, cross | 8 / 123 / 0.00% | 8 / 0.125 / 123 |
| 64/16, $10^4$, cross | 21 / 65 / 20.37% | 18 / 0.281 / 54 |
| 128/16, $10^4$, cross | 43 / 242 / 4.76% | 38 / 0.297 / 231 |
| 128/16, $10^2$, cross | 32 / 293 / 29.65% | 26 / 0.203 / 226 |
| 128/16, $10^6$, cross | 64 / 273 / 3.80% | 62 / 0.484 / 263 |
| 128/16, $10^4$, ring | 43 / 237 / 3.95% | 52 / 0.406 / 228 |

七个设计问题上，adaptive 的平均/最大 gap 为 8.93%/29.65%。三个参数冻结后的验证
问题上，平均/最大 gap 为 27.34%/54.20%；88/11、$10^4$ diagonal，120/15、$10^2$
branching，以及 152/19、$10^6$ parallel 上 adaptive/oracle 分别为 278/243、
228/201 和 367/238 次。设计组与验证组分开报告，以保持留出检验的独立性。

## 9. Seed、RHS 与中心问题计时

五个系数 seed 和同一矩阵上的六种 RHS 汇总为：

| 试验 | 方法 | 收敛/慢/发散 | 循环总和 | 平均 | 最坏 |
|---|---|---:|---:|---:|---:|
| 5 seed | adaptive | 5/0/0 | 327 | 65.40 | 75 |
| 5 seed | global-reference | 5/0/0 | 866 | 173.20 | 321 |
| 6 RHS | adaptive | 6/0/0 | 469 | 78.17 | 90 |
| 6 RHS | global-reference | 6/0/0 | 1507 | 251.17 | 328 |
| 6 RHS | geometric | 6/0/0 | 156967 | 26161.17 | 29310 |

中心 128/16、$10^4$ cross-channel 问题上的五次预热后重复计时为：

| 方法 | setup 中位数 ms（Q1--Q3） | solve 中位数 ms（Q1--Q3） | total 中位数 ms（Q1--Q3） | 循环 |
|---|---:|---:|---:|---:|
| adaptive | 375.720（371.005--379.193） | 225.587（214.829--227.634） | 598.639（590.434--601.306） | 242 |
| global-reference | 1903.654（1888.064--1931.348） | 7252.850（7194.972--7316.201） | 9156.505（9083.036--9256.025） | 3227 |

adaptive 在中心问题选择 $m=43$；global-reference 的 setup、循环数和完整求解时间
均明显更高，total 中位数之比为 15.30。

## 10. 复现

在 `code` 目录执行：

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

脚本优先使用 CMake，并提供等价的严格 C++17 直接构建路径。四份完整逐行结果位于
`code/results`；`full` 生成正式结果，`quick` 用于开发检查并默认写入
`code/build/quick-results`，不会覆盖正式结果。
