# two_grids_iteration v5.7.0

## 核心模块

| 模块 | 内容 |
|---|---|
| `energy_interpolation.hpp` | 几何插值与全局能量数值参考 |
| `global_pcg.hpp` | 可延续的全局 PCG 路径与尺度自适应 fast/reuse 选择 |
| `two_grid_solver.hpp` | Galerkin 粗算子、对称平滑、收敛状态与收敛因子 |

实现采用固定粗点、精确 Galerkin 粗解和对称 Gauss--Seidel 的两网格结构。

## 尺度自适应策略

令 `n=grid.intervals()`、`n_H=n/grid.ratio()`。fast 面向单个或少量 RHS，只构造一个
候选：`n_H<=8` 时取 `n/8`；其余情况按矩阵对角尺度比 `<1e3`、`[1e3,1e5)`、
`>=1e5` 分别取 `n/4`、`n/3`、`n/2`。它不运行 pilot，setup 只有一次 PCG 路径推进、
插值组装和粗算子构造。

reuse 面向多 RHS，取 `n/8,3n/16,n/4,n/3,n/2` 五个归一化候选。每个候选最多运行
`n/2` 个 pilot 循环，以后半段
残量的几何平均收缩率预测完整循环数；非有限或非收缩的轨迹按正式循环上限处理。
候选按路径顺序流式处理，只保留当前优胜者。

离线 oracle 额外报告所选步数的归一化位置 `m/(1/h)`，用于检验理论中的尺度窗口；
该量不进入在线选择。

```cpp
tgi::AdaptiveGlobalPcgOptions options;
options.policy = tgi::AdaptiveGlobalPcgPolicy::Fast;
// options.policy = tgi::AdaptiveGlobalPcgPolicy::Reuse;
```

## 构建与验证

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

脚本优先使用 CMake；没有 CMake 时采用等价的严格 C++17 直接构建。可通过
`TGI_THREADS`、`TGI_BUILD_DIR` 和 `TGI_STEP_TIMEOUT_SECONDS` 调整运行。

## 实验

| 入口 | 内容 |
|---|---|
| `experiment1_two_grid_comparison` | 尺寸、对比度、六类拓扑和 256/16 大尺度比较 |
| `experiment2_finite_path` | 能量单调下降与循环数非单调变化的直接证据 |
| `experiment3_oracle_validation` | 七个设计问题和三个冻结后验证问题的归一化 step--2 离线 oracle |
| `experiment4_submission_robustness` | 五种 seed、六类 RHS 与中心问题重复计时 |

正式求解从零初值运行到相对残量 `1e-6`，循环上限为 20000。结果同时报告
`converged`、`slow-limit`、`diverged`、全程有效收敛因子和末端收敛因子。
中心计时表仅在三种被比较方法均达到正式容差后生成。
