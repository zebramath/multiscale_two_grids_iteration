# two_grids_iteration v7.3.0

## 核心模块

| 模块 | 内容 |
|---|---|
| `energy_interpolation.hpp` | 几何插值与全局能量数值参考 |
| `global_pcg.hpp` | 可延续的全局 PCG 路径与尺度自适应有限步构造 |
| `two_grid_solver.hpp` | Galerkin 粗算子、对称平滑、统一稳态迭代状态与收敛因子 |
| `multilevel_solver.hpp` | 多层层次校验、对称递归 V-cycle、算子与插值复杂度 |

实现采用固定粗点、精确 Galerkin 粗解和对称 Gauss--Seidel 的两网格结构。

## 理论知情的轻量经验策略

adaptive 的定位为

$$
\boxed{\text{theory-informed lightweight empirical selector}}.
$$

令 `n=grid.intervals()`、`n_H=n/grid.ratio()`。adaptive 只构造一个
候选：`n_H<=8` 时取 `n/8`；其余情况按矩阵对角尺度比 `<1e3`、`[1e3,1e5)`、
`>=1e5` 分别取 `n/4`、`n/3`、`n/2`。setup 只有一次 PCG 路径推进、
插值组装和粗算子构造。

理论只支持固定能量衰减的 `m=O(1/h)` 充分尺度。`1/8,1/4,1/3,1/2` 与
`1e3,1e5` 均为当前受控问题族上的经验参数，不宣称最优或 near-oracle。离线 oracle
仅在预设 step--2 窗口内作描述性比较；它不进入在线选择，也不提供参数最优性保证。

```cpp
auto result = tgi::build_adaptive_global_pcg_interpolation(
    grid, matrix, geometric.prolongation, 4);
```

最后一个参数只控制并行线程数；在线策略固定为一次路径推进和单候选构造。

## 构建与验证

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh supplemental
./scripts/run_validation.sh multilevel
./scripts/run_validation.sh full
```

脚本优先使用 CMake；没有 CMake 时采用等价的严格 C++17 直接构建。可通过
`TGI_THREADS`、`TGI_BUILD_DIR` 和 `TGI_STEP_TIMEOUT_SECONDS` 调整运行。`supplemental`
只运行实验 5--6；`multilevel` 只运行实验 7。`quick`
结果默认写入构建目录下的 `quick-results`，不会覆盖 `results` 中的正式完整结果；
可分别通过 `TGI_QUICK_RESULTS_DIR` 和 `TGI_RESULTS_DIR` 改写输出位置。

## 实验

| 入口 | 内容 |
|---|---|
| `experiment1_two_grid_comparison` | 尺寸、对比度、六类拓扑和 256/16 大尺度比较 |
| `experiment2_finite_path` | 能量单调下降与循环数非单调变化的直接证据 |
| `experiment3_oracle_validation` | 七个设计问题和三个冻结后验证问题的归一化 step--2 窗口受限离线 oracle |
| `experiment4_submission_robustness` | 五种 seed、六类 RHS 与中心问题重复计时 |
| `experiment5_stopping_ablation` | adaptive 与固定归一化步数、逐列固定残差停止的消融 |
| `experiment6_fixed_physical_refinement` | 固定物理系数场的三层嵌套网格加密 |
| `experiment7_multilevel_pilot` | 三层 V-cycle 与首层精确两网格的配对比较 |

正式求解从零初值运行到相对残量 `1e-6`；adaptive/reference 的循环上限为 20000，
geometric 为 30000。结果同时报告
`converged`、`slow-limit`、`diverged`、全程有效收敛因子和末端收敛因子。
中心计时表仅在两种被比较方法均达到正式容差后生成。

实验 5 的三种策略从同一几何插值出发并求解同一组全局 PCG 列方程；结果以确定性的
列迭代总数和两网格循环数比较，不用单次墙钟时间作结论。实验 6 将通道物理宽度固定为
`1/16`、背景划分固定为 `8 x 8`，并在求解前检查相邻嵌套网格所有共享节点的系数值。
实验 7 在两个 Galerkin 层次上独立构造 transfer，使用一次前向/后向 Gauss--Seidel
和最粗层精确 Cholesky；它报告确定性的循环数、算子复杂度 $C_A$ 与插值复杂度 $C_P$，
不用单次墙钟时间作结论。
