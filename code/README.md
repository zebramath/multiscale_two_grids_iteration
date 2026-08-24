# two_grids_iteration v4.1.0

## 核心模块

| 模块 | 内容 |
|---|---|
| `energy_interpolation.hpp` | 几何、局部和全局能量插值 |
| `global_pcg.hpp` | 增量全局 PCG 路径和渐进式候选竞速 |
| `two_grid_solver.hpp` | Galerkin 粗算子、对称平滑和两网格循环 |
| `support_expansion.hpp` | 局部支撑研究基础设施 |

v4.1 不包含递归多层层次、V-cycle 或多层预条件 CG。

## 自适应算法

若几何候选不够好，选择器沿同一增量 PCG 路径广筛
`m=12,20,28,36,44,52,60`，每个候选运行至多 64 次两网格 pilot。随后确认预测最好
的三个候选、最小步锚点和 3% 松弛下的稀疏优先候选，确认上限为 384 次；最后检查
优胜者的 step--2 邻点和它与次优者之间的中点。

候选规则不读取拓扑、对比度、种子或 oracle 步数。选择墙钟时间包括所有候选两网格
层次的建立和 pilot；主实验的 Setup 还包括几何初值构造。

## 构建与验证

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

可通过 `TGI_THREADS` 和 `TGI_STEP_TIMEOUT_SECONDS` 调整线程数与单步超时。

## 三个实验

| 入口 | 内容 |
|---|---|
| `experiment1_two_grid_comparison` | 单因素尺寸、对比度、六拓扑两网格主比较 |
| `experiment2_finite_path` | 三个 cross-channel 问题的有限路径机制证据 |
| `experiment3_oracle_validation` | 七问题 step--2 离线 oracle 选择质量 |

每个入口生成一份同名 TXT，不生成 CSV。所有公开循环数均由候选独立求解到相对残差
`1e-6` 得到。
