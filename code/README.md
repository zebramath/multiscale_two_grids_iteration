# two_grids_iteration v4.7.0

## 核心模块

| 模块 | 内容 |
|---|---|
| `energy_interpolation.hpp` | 几何插值与全局能量数值参考 |
| `global_pcg.hpp` | 增量全局 PCG 路径和 v4.7 单参数预算选择器 |
| `two_grid_solver.hpp` | Galerkin 粗算子、对称平滑和两网格循环 |

代码采用固定粗点的对称两网格结构，核心路径由全局/有限能量插值、Galerkin 粗校正和
两档预算选择组成。实验入口直接承载对应研究逻辑，模块边界与论文中的算法对象一致。

## 自适应策略

以代表性右端项的短两网格残差尾预测完整循环数。公开参数
`expected_rhs_count=R` 连续控制候选间距、最大路径步数、pilot 上限和近优容忍度；
较大 `R` 还允许在预测优胜者的 `m±2` 处做至多两次细化。默认 `R=1` 只评估
`m=0,12,32,52` 且每个候选的 pilot 上限为 16；`R=256` 评估
`m=0,12,20,...,60` 且上限为 160。候选达到正式容差后提前停止，尺度、对比度和拓扑
变化通过候选的实际残量历史影响选择。

典型设置：

```cpp
tgi::AdaptiveGlobalPcgOptions options;
options.expected_rhs_count = 1.0;    // 低 setup
// options.expected_rhs_count = 256; // 多右端项，偏向低循环数
```

## 构建与验证

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

脚本优先使用 CMake；没有 CMake 时自动使用 C++17 严格直接构建。可通过
`TGI_THREADS`、`TGI_BUILD_DIR` 和 `TGI_STEP_TIMEOUT_SECONDS` 调整运行。

## 四个实验

| 入口 | 内容 |
|---|---|
| `experiment1_two_grid_comparison` | 尺寸、对比度、六拓扑的两网格主比较 |
| `experiment2_finite_path` | 三个问题上的能量单调/循环非单调证据 |
| `experiment3_oracle_validation` | 七问题 step--2 离线 oracle 质量 |
| `experiment4_submission_robustness` | 五种系数种子、六类 RHS 迁移和五次重复计时 |

每个入口生成同名 TXT。正式循环数均从零初值独立求解到相对残差 `1e-6`，循环上限为
12000。输出区分 `converged`、`slow-limit` 与 `diverged`；计时汇总采用共同收敛子集。
