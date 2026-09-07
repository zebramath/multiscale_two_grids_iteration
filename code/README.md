# two_grids_iteration v8.4.0

## 实现

| 模块 | 内容 |
|---|---|
| `linear_algebra.hpp` | CSR 稀疏矩阵、并行乘法和稀疏 Cholesky |
| `diffusion_problem.hpp` | 规则网格、高对比系数场和扩散算子 |
| `energy_interpolation.hpp` | 几何插值与能量数值参考 |
| `global_pcg.hpp` | 可延续的全局 PCG 路径与自适应有限步规则 |
| `two_grid_solver.hpp` | Galerkin 粗算子、对称两网格循环和收敛指标 |
| `multilevel_solver.hpp` | 递归 V-cycle、算子复杂度和插值复杂度 |

实现采用固定粗点、精确 Galerkin 粗解和一次前向/后向 Gauss--Seidel。`adaptive` 根据
粗网格分辨率与矩阵对角尺度比，在 `n/8`、`n/4`、`n/3` 和 `n/2` 中选择一个 PCG
路径位置。标准 PCG 估计在固定对比度下给出固定比例误差压缩的 `O(1/h)` 充分尺度；
进入理论中能量极小值点的局部区域还需要对数因子。具体比例与阈值由设计问题组确定，
不宣称最优。

```cpp
auto initial = tgi::build_geometric_interpolation(grid);
auto result = tgi::build_adaptive_global_pcg_interpolation(
    grid, matrix, initial.prolongation, 4);
```

## 构建与运行

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh supplemental
./scripts/run_validation.sh multilevel
./scripts/run_validation.sh endpoint
./scripts/run_validation.sh full
```

脚本优先使用 CMake；无 CMake 时采用严格 C++17 直接构建；完整模式绘图需要 Python 3
与 Matplotlib。环境变量 `TGI_THREADS`、
`TGI_BUILD_DIR`、`TGI_RESULTS_DIR`、`TGI_QUICK_RESULTS_DIR` 和
`TGI_STEP_TIMEOUT_SECONDS` 可控制线程、目录与常规步骤超时。`quick` 构建全部入口并运行
实验 1 的三个代表问题；`supplemental` 运行实验 5--6；`multilevel` 运行实验 7；
`endpoint` 仅运行实验 8，并让两种方法持续迭代至收敛而不设循环数上限或脚本超时；
`full` 运行八组正式实验。

## 四个数值主题与八个实验入口

| 主题 | 入口 | 内容 |
|---|---|---|
| 两网格性能与逐步扫描 | `experiment1_two_grid_comparison` | 尺度、对比度、六类拓扑和 256/16 扩展问题 |
| 两网格性能与逐步扫描 | `experiment2_step_scan` | cross 与 winding-ring 上逐一扫描 `m=1,...,128`，输出归一化能量与 $\rho_{\mathrm{eff}}$ |
| 两网格性能与逐步扫描 | `experiment8_interpolation_endpoint` | 中心问题上几何插值与能量极小插值的无循环上限单次对照 |
| 规则质量、稳健性与代价 | `experiment3_oracle_validation` | 设计组与冻结验证组上的窗口受限离线采样参考 |
| 规则质量、稳健性与代价 | `experiment4_submission_robustness` | 五个 seed、六类 RHS 与中心问题重复计时 |
| 规则质量、稳健性与代价 | `experiment5_stopping_ablation` | 自适应有限步、固定步数与固定列残量消融 |
| 固定物理系数场加密 | `experiment6_fixed_physical_refinement` | 固定物理系数场的三层嵌套加密 |
| 多层初步验证 | `experiment7_multilevel_pilot` | 三层 V-cycle 与首层精确两网格配对 |

正式求解从零初值开始，相对欧氏残量容差为 `1e-6`。常规循环上限为 20000；Exp2
逐步扫描使用 12000 个循环的固定观测上限；未在该上限前收敛的点以第 12000 次循环的
残量计算经验因子。常规结果报告循环数、最终残量、经验有效收敛因子
$\rho_{\mathrm{eff}}$、末端因子、插值密度和收敛状态。Exp2 同时记录
$J(W)=\tfrac12\operatorname{tr}(P^\top A P)$ 与 $\rho_{\mathrm{eff}}$；前者沿路径单调，
后者是实际残量历程指标，不等同于理论能量范数因子 $\rho_{TG}$。Exp4 的计时比较由几何
插值出发，双方采用同一 Jacobi--PCG 路径，报告预热后五次测量的算术平均时间。

实验 5 以列迭代总数衡量确定性 setup 工作量。实验 6 固定物理通道宽度和背景分区，并
逐点核对嵌套网格共享节点。实验 7 在两个 Galerkin 转移上独立构造插值，报告
$C_A$、$C_P$ 及 V-cycle/两网格循环数比。实验 8 固定中心 128/16、对比度 $10^4$、
cross-channel 和常数右端项，仅更换插值矩阵；几何插值与列残量容差 $10^{-10}$ 的能量
极小插值分别需要 85524 和 3227 个循环达到 $10^{-6}$，循环数比为 26.503。
