# 多尺度有限 Krylov 能量插值研究方案 v5.4

## 1. 研究定位

本项目研究高对比扩散问题中，全局能量极小插值的有限 PCG 路径与实际两网格收敛
目标之间的非单调错位。中心现象是：插值能量随 PCG 迭代持续下降，而固定磨光器下的
两网格循环数可以先快速下降、在有限步达到最小、随后回升，并在高精度能量终点显著
增大。

研究集中回答三个问题：

1. 能量距离、粗空间几何和两网格谱因子之间有什么严格联系；
2. 有限 Krylov 误差在什么条件下改善磨光后最坏误差模态；
3. 最优有限步在什么条件下具有 $\Theta(h^{-1})$ 尺度，以及如何据此在少 RHS 与多 RHS
   场景间平衡 setup 和求解成本。

二维规则网格、固定粗点、对称 Gauss--Seidel、精确 Galerkin 粗校正和两网格结构构成
受控研究平台，使有限 Krylov 插值对粗空间质量的作用能够被单独识别。

### 1.1 相关工作与本文差异

能量极小插值从能量、有约束近核空间再现和稀疏支撑构造延伸到并行实现，形成了成熟的
AMG 插值路线（[Wan--Chan--Smith 2000](https://doi.org/10.1137/S1064827598334277)、
[Olson--Schroder--Tuminaro 2011](https://doi.org/10.1137/100803031)、
[Janna 等 2023](https://doi.org/10.1137/22M1513794)）。另一方面，直接以两网格收敛率
定义最优插值说明能量代理与求解目标应明确区分
（[Brannick 等 2018](https://doi.org/10.1137/17M1123456)）。参数化数值方法的近期研究
则采用“真实目标、廉价代理、选择开销、独立验证”的闭环
（[Frasca-Caccia--Singh 2023](https://doi.org/10.1007/s10915-023-02324-0)）。

本项目沿用这一研究范式，但研究对象是同一能量方程的有限 Krylov 路径：先证明能量
终点与两网格目标的结构性错位，再以尺度化候选和短两网格轨迹选择有限步，并分别报告
候选覆盖误差、oracle gap、setup--solve 成本和跨 RHS 迁移。

## 2. 模型问题与比较对象

考虑

$$
-\nabla\cdot(a(x)\nabla u)=f\quad\text{in }\Omega=(0,1)^2,
\qquad u=0\quad\text{on }\partial\Omega,
$$

其中 $a(x)>0$，对比度取 $10^2$、$10^4$ 或 $10^6$。节点型有限差分离散得到
对称正定系统 $A_hu_h=b_h$。系数结构包括交叉、曲折、对角、平行、分支和闭合曲环
六类高导通道，并叠加由 seed 控制的块状高导背景。

按细点和粗点分块，写成

$$
A=\begin{bmatrix}B&C\\ C^T&A_{CC}\end{bmatrix},
\qquad
P(W)=\begin{bmatrix}W\\I\end{bmatrix}.
$$

比较四种方法：

1. `adaptive-fast`：尺度自适应的单候选有限 PCG 插值；
2. `adaptive-reuse`：在归一化路径区间内由短两网格试验选择有限 PCG 插值；
3. `global-reference`：全局能量方程达到相对欧氏残量 $10^{-10}$ 的高精度数值参考；
4. `geometric`：规则粗点上的双线性几何插值 $P_G$。

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

PCG 路径 $W_m$ 在 $B$ 能量范数下单调逼近 $W_*$。两网格误差传播算子

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

因此有限步改善条件不依赖一对任意选取的单奇异向量，而要求 Krylov 方向在整个主导
子空间上产生正定对齐；标量公式是重数为一时的特例。

全局稳定界还能把固定性能优势变成定量能量屏障。若候选相对能量终点的
$\sqrt{\rho_{TG}}$ 改善至少为 $\gamma$，则

$$
J(W)-J(W_*)\ge
\frac{\lambda_{\min}(S)}2
\frac{\gamma^2}{\lVert T\rVert_2^2-\gamma^2}.
$$

因此具有固定幅度优势的有限步候选不能无限逼近能量终点，这把非单调现象从方向性条件
加强为可分辨的能量距离结论。

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

这给出 $m_{\mathrm{opt}}(h)=\Theta(h^{-1})$ 的条件性估计。上界直接来自标准 PCG
条件数界，不再另设路径指数上包络；额外路径条件只负责必要的下界。常数允许随固定的
对比度范围、拓扑类别和粗化比变化；该结论描述两网格性能最优点，而 PCG 的代数终止
上界仍由细点系统维数决定。七个离线 step--2 sampled oracle 的
$m_{\mathrm{or}}/(1/h)$ 位于
$0.125$--$0.484$，支持当前 $[1/8,1/2]$ 归一化窗口；这些离散采样点是尺度证据，
不等同于连续候选路径上的 $m_{\mathrm{opt}}(h)$。

若归一化路径投影具有与网格无关的连续模 $\omega_\Pi$，且
$\lVert T_h\rVert_2\le M_T$，记舍入后候选的真实覆盖半径为
$\delta_{\mathcal S,h}\le\delta_\mathcal S+h/2$，则令
$\varepsilon_h=M_T\omega_\Pi(\delta_{\mathcal S,h})$ 后有

$$
\sqrt{\rho_{TG}(\widehat m_h)}-\sqrt{\rho_{TG}(m_{\mathrm{opt}})}
\le\varepsilon_h,
$$

以及更紧的平方目标界

$$
\rho_{TG}(\widehat m_h)-\rho_{TG}(m_{\mathrm{opt}})
\le\varepsilon_h\min\left\{2M_T,
2\sqrt{\rho_{TG}(m_{\mathrm{opt}})}+\varepsilon_h\right\}.
$$

该结论无需 Lipschitz 正则性，也无需舍入点仍位于原区间。固定候选数的离散损失由
归一化覆盖半径而不是网格规模控制。当前 reuse 集合在 $[1/8,1/2]$ 上的连续覆盖
半径为 $1/12$。

### 4.2 Fast：低 setup 单候选

令 $n_H=1/H$。fast 直接取

$$
m_{\mathrm{fast}}=
\begin{cases}
\operatorname{round}(n/8),&n_H\le8,\\
\operatorname{round}(n/4),&n_H>8,\ \chi_A<10^3,\\
\operatorname{round}(n/3),&n_H>8,\ 10^3\le\chi_A<10^5,\\
\operatorname{round}(n/2),&n_H>8,\ \chi_A\ge10^5.
\end{cases}
$$

该档只推进一次 PCG 路径、组装一次插值并建立一次粗算子，不运行 pilot。规则只使用
粗空间分辨率和矩阵对角尺度比；$10^3$ 与 $10^5$ 是三个测试对比度的对数中点。

### 4.3 Reuse：五个归一化候选

reuse 使用

$$
\mathcal M_{\mathrm{reuse}}
=\operatorname{round}\left\{n/8,3n/16,n/4,n/3,n/2\right\},
\qquad L=\operatorname{round}(n/2).
$$

候选数量固定为五个，搜索区间和 pilot 长度均随分辨率伸缩。几何端点退出在线选择，
$3n/16$ 细化 oracle 较集中的左半区间。实现沿同一条增量 PCG 路径依次生成候选，并
流式保留当前优胜者。

对候选 $m$，记 pilot 实际运行 $\ell_m\le L$ 次，第 $k$ 次后的相对残量为 $r_{m,k}$，
$k_0=\lfloor\ell_m/2\rfloor$。估计

$$
\widehat\rho_m=
\left(\frac{r_{m,\ell_m}}{r_{m,k_0}}\right)^{1/(\ell_m-k_0)},
$$

并外推达到正式容差所需的循环数 $\widehat N_m$。reuse 选择预测循环数最小者，预测
相同时取较小 $m$。令 $e_m=\widehat N_m-N_m$，$m_{\mathcal M}$ 为候选集内真实
循环数最少者，则无需控制其他候选便有精确的后悔上界

$$
0\le N_{\widehat m}-N_{m_{\mathcal M}}
\le e_{m_{\mathcal M}}-e_{\widehat m}.
$$

若真实最优候选至多被高估 $\varepsilon_+^{\mathrm{opt}}$，实际所选候选至多被低估
$\varepsilon_-^{\mathrm{sel}}$，则加性损失至多为二者之和。相应的单侧相对误差给出

$$
\frac{N_{\widehat m}}{N_{m_{\mathcal M}}}
\le\frac{1+\delta_+^{\mathrm{opt}}}
{1-\delta_-^{\mathrm{sel}}}.
$$

统一双侧误差界只是上述结论的特例；新形式减弱了条件，并准确指出 pilot 误差中真正
决定选择损失的两部分。

### 4.4 多 RHS 成本

若策略 $a,b$ 的一次 setup 为 $S_a,S_b$，单 RHS 求解时间为 $T_a,T_b$，则

$$
C_a(R)=S_a+RT_a,
\qquad
R_{\mathrm{break}}=\frac{S_b-S_a}{T_a-T_b}.
$$

当 $S_b>S_a$ 且 $T_b<T_a$ 时，$R>R_{\mathrm{break}}$ 才使策略 $b$ 的总成本更低。
fast 面向低 setup 的少 RHS 工作量，reuse 面向能够摊销选择成本的多 RHS 工作量。

## 5. 实验设计与指标

所有方法使用相同矩阵、粗点、右端项、前后 Gauss--Seidel 和 Galerkin 粗算子。正式
求解从零初值开始，相对欧氏残量容差为 $10^{-6}$，循环上限为 20000。

核心证据分为四组：

1. 13 个问题比较尺寸、对比度和通道拓扑，其中两组 256/16 作为大尺度扩展；
2. 3 条归一化有限 PCG 路径比较能量超额与两网格循环数；
3. 7 个设计问题和 3 个参数冻结后验证问题上的 step--2 离线 oracle；
4. 5 个系数 seed、6 个迁移 RHS 和中心 128/16 问题的 5 次重复计时。

主要指标包括循环数、收敛状态、插值密度，以及以实际执行循环数 $k$ 定义的全程有效
收敛因子

$$
\rho_{\mathrm{eff}}=r_k^{1/k},
$$

以及最后 32 次循环的尾部因子。达到上限且残量仍收缩的结果记为 `slow-limit`；残量
非有限或持续增长的结果记为 `diverged`。时间只在中心 128/16 问题上比较三个完整
收敛的方法，并报告预热后五次测量的 Q1、中位数和 Q3。

## 6. 跨尺寸、对比度与拓扑的主结果

| 轴 | $1/h,1/H$ | 对比度/拓扑 | fast `m/cycles` | reuse `m/cycles` | reference cycles | geometric |
|---|---|---|---:|---:|---:|---:|
| size | 32,8 | $10^4$/cross | 4 / 37 | 4 / 37 | 54 | 40 |
| size | 64,8 | $10^4$/cross | 8 / 123 | 8 / 123 | 264 | 302 |
| size | 64,16 | $10^4$/cross | 21 / 65 | 16 / 60 | 321 | 20000 slow |
| center | 128,16 | $10^4$/cross | 43 / 242 | 43 / 242 | 3227 | 20000 slow |
| large | 256,16 | $10^4$/cross | 85 / 998 | 85 / 998 | 13948 | 20000 slow |
| large | 256,16 | $10^4$/ring | 85 / 1002 | 128 / 950 | 3950 | 20000 slow |
| contrast | 128,16 | $10^2$/cross | 32 / 293 | 24 / 237 | 756 | 1501 |
| contrast | 128,16 | $10^6$/cross | 64 / 273 | 64 / 273 | 3471 | 20000 slow |
| topology | 128,16 | $10^4$/meandering | 43 / 266 | 43 / 266 | 759 | 20000 slow |
| topology | 128,16 | $10^4$/diagonal | 43 / 346 | 43 / 346 | 418 | 20000 slow |
| topology | 128,16 | $10^4$/parallel | 43 / 260 | 43 / 260 | 421 | 20000 slow |
| topology | 128,16 | $10^4$/branching | 43 / 317 | 43 / 317 | 615 | 20000 slow |
| topology | 128,16 | $10^4$/ring | 43 / 237 | 43 / 237 | 768 | 20000 slow |

fast、reuse 和 global-reference 均在 13 个问题上收敛，累计循环数分别为 4459、4346
和 28972；global-reference 分别为两种有限策略的 6.50 倍和 6.67 倍。geometric 在
3 个问题上收敛，其余 10 个为 `slow-limit`，没有发散结果。

尺度伸缩在 256/16 上保持有效：cross-channel 的 fast/reuse 均为 998 次，而 reference
为 13948 次；winding-ring 的 reuse 为 950 次，reference 为 3950 次。中心问题的有效
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

## 8. 离线 sampled oracle 选择质量

离线 sampled oracle 检查 $m=0$ 以及从 $n/8$ 到 $n/2$ 的 step--2 路径点。七个设计
问题用于确定归一化区间，三个参数冻结后的新尺度/拓扑/对比度组合用于独立检查。超过
6000 次仍未达到容差的候选退出最优比较，因为它们已经无法优于表中的有限最小值。
下表的 oracle 表示该离散评估集合中的最优点。

| 问题 | fast `m/cycles/gap` | reuse `m/cycles/gap` | oracle `m/(1/h)/cycles` |
|---|---:|---:|---:|
| 32/8, $10^4$, cross | 4 / 37 / 0.00% | 4 / 37 / 0.00% | 4 / 0.125 / 37 |
| 64/8, $10^4$, cross | 8 / 123 / 0.00% | 8 / 123 / 0.00% | 8 / 0.125 / 123 |
| 64/16, $10^4$, cross | 21 / 65 / 20.37% | 16 / 60 / 11.11% | 18 / 0.281 / 54 |
| 128/16, $10^4$, cross | 43 / 242 / 4.76% | 43 / 242 / 4.76% | 38 / 0.297 / 231 |
| 128/16, $10^2$, cross | 32 / 293 / 29.65% | 24 / 237 / 4.87% | 26 / 0.203 / 226 |
| 128/16, $10^6$, cross | 64 / 273 / 3.80% | 64 / 273 / 3.80% | 62 / 0.484 / 263 |
| 128/16, $10^4$, ring | 43 / 237 / 3.95% | 43 / 237 / 3.95% | 52 / 0.406 / 228 |

七个设计问题上，fast 的平均/最大 gap 为 8.93%/29.65%，对应一次候选构造和零
pilot；reuse 为 4.07%/11.11%，对应固定五个候选。三个参数冻结后的验证问题上，fast
的平均/最大 gap 为 27.34%/54.20%，reuse 为 19.44%/43.62%。具体结果为：88/11、
$10^4$ diagonal 上 fast/reuse/oracle 为 278/349/243 次，120/15、$10^2$ branching
上为 228/228/201 次，152/19、$10^6$ parallel 上为 367/241/238 次。reuse 在后两例
达到或接近 fast，并在最大尺度验证例接近 oracle；第一例的短轨迹外推选到右端点，说明
有限 pilot 预测仍是策略的主要误差来源。设计组与验证组分开报告，避免以留出结果继续
调整规则。

## 9. Seed、RHS 与中心问题计时

五个系数 realization 的汇总为：

| 方法 | 收敛/慢/发散 | 循环总和 | 平均 | 最坏 |
|---|---:|---:|---:|---:|
| adaptive-fast | 5/0/0 | 327 | 65.40 | 75 |
| adaptive-reuse | 5/0/0 | 278 | 55.60 | 60 |
| global-reference | 5/0/0 | 866 | 173.20 | 321 |

候选只由常数 RHS 构造一次，再用于六种 RHS。汇总为：

| 方法 | 收敛/慢/发散 | 循环总和 | 平均 | 最坏 |
|---|---:|---:|---:|---:|
| adaptive-fast | 6/0/0 | 469 | 78.17 | 90 |
| adaptive-reuse | 6/0/0 | 336 | 56.00 | 60 |
| global-reference | 6/0/0 | 1507 | 251.17 | 328 |
| geometric | 1/5/0 | 19120 | 19120.00 | 19120 |

中心 128/16、$10^4$ cross-channel 问题上的五次预热后重复计时为：

| 方法 | setup 中位数 ms（Q1--Q3） | solve 中位数 ms（Q1--Q3） | total 中位数 ms（Q1--Q3） | 循环 |
|---|---:|---:|---:|---:|
| adaptive-fast | 349.049（342.178--355.661） | 211.420（198.172--216.587） | 552.731（547.221--562.297） | 242 |
| adaptive-reuse | 846.899（844.134--869.044） | 209.807（205.318--212.680） | 1056.813（1056.706--1071.007） | 242 |
| global-reference | 1819.059（1759.163--1839.425） | 6770.640（6752.234--6907.390） | 8608.713（8589.699--8666.554） | 3227 |

fast 与 reuse 在中心问题选出相同的 $m=43$ 插值，两者求解时间差属于重复计时波动；
fast 省去候选 pilot，因此具有更低的 setup 和端到端时间。global-reference 的 setup、
循环数和完整求解时间均明显更高。

## 10. 复现

在 `code` 目录执行：

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

脚本优先使用 CMake，并提供等价的严格 C++17 直接构建路径。四份完整逐行结果位于
`code/results`；`full` 生成正式结果，`quick` 用于开发检查。
