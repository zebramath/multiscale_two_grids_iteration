# two_grids_iteration v4.2.0

## 核心模块

| 模块 | 内容 |
|---|---|
| `energy_interpolation.hpp` | 几何、局部和全局能量插值 |
| `global_pcg.hpp` | 增量全局 PCG 路径和 v4.2 可调选择器 |
| `two_grid_solver.hpp` | Galerkin 粗算子、对称平滑和两网格循环 |

代码不包含递归多层层次、V-cycle 或多层预条件 CG。未进入研究主线的支撑扩张模块
和旧版多阶段候选路由已删除。

## 自适应策略

默认候选为 `m=0,12,20,...,60`，以代表性右端项的短两网格残差尾预测完整循环数。
公开参数 `expected_rhs_count=R` 连续控制 pilot 长度和近优容忍度；较大 `R` 还允许在
预测优胜者的 `m±2` 处做至多两次细化。策略不读取尺寸、对比度、拓扑或 oracle。

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

## 三个实验

| 入口 | 内容 |
|---|---|
| `experiment1_two_grid_comparison` | 尺寸、对比度、六拓扑的两网格主比较 |
| `experiment2_finite_path` | 三个问题上的能量单调/循环非单调证据 |
| `experiment3_oracle_validation` | 七问题 step--2 离线 oracle 质量 |

每个入口生成同名 TXT。正式循环数均从零初值独立求解到相对残差 `1e-6`，不复用
候选 pilot 状态。
