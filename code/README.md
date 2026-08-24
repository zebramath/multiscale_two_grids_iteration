# two_grids_iteration v4.0.0

## 目录结构

`multigrid` 包含六个核心模块：

| 模块 | 内容 |
|---|---|
| `energy_interpolation.hpp` | 几何、局部、全局能量插值，强连接距离插值和参考剪枝 |
| `support_expansion.hpp` | 强连接支撑扩张和残差预算自适应支撑 |
| `global_pcg.hpp` | 增量全局 PCG 路径和自适应检查点选择 |
| `two_grid_solver.hpp` | Galerkin 粗算子、对称平滑、粗解和两网格循环 |
| `multilevel_solver.hpp` | 递归对称 V-cycle、独立多层迭代和预条件 CG |
| `multilevel_hierarchy.hpp` | 几何、自适应 PCG、全局能量极小三种多层构造策略 |

`CoarseSetupReport` 保留能量和层次统计供内部一致性检查。v4.0 的公开主实验不再
单列复杂度/谱诊断扫描，只在有限路径实验中报告验证中心现象所需的能量超额、密度
和实际循环数。

`experiments` 只包含四个可执行入口，每个入口生成一份完整 TXT。

## 自适应算法

选择器复用一条增量式全局能量 PCG 路径，最多建立四个层次：几何初值、最小步锚点、区间中点和一个补充点。每个候选运行至多 24 个两网格 pilot，不在 setup 阶段求解到完整容差。候选位置不读取系数场名称、对比度、种子或预设最优步数。

几何初值或最小步锚点已经足够好时立即停止；否则检查中点。中点仍较难时，使用归一化 PCG 能量残差决定检查锚点—中点内点或最大预算端点。在最小预测的 13% 范围内选择最小 `m`。v4.0 用最大预算端点替代原来的上半区中点，避免强非渐近问题在有效区间之前停止。

多层构造在首层使用实验指定的 $H$，此后按因子 2 递归粗化至单个未知量。三种策略在每层保持一致：几何插值、低预算自适应有限 PCG、或全局能量极小插值。所得对称 V-cycle 既可以直接定常迭代，也可以作为标准 PCG 预条件子。

## 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

## 验证

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

可通过 `TGI_THREADS` 和 `TGI_STEP_TIMEOUT_SECONDS` 调整线程数与单步超时。

## 四个整合实验

| 入口 | 内部测试角度 |
|---|---|
| `experiment1_two_grid_comparison` | 12 个尺寸/对比度/拓扑问题上的自适应、global-exact、几何两网格比较 |
| `experiment2_finite_path` | 3 个代表问题上能量单调而两网格性能非单调的紧凑证据 |
| `experiment3_oracle_validation` | 8 问题 step-2 离线 oracle 选择质量 |
| `experiment4_multilevel_comparison` | 12 个问题上的多层 V-cycle 与多层预条件 CG 比较 |

每个入口只生成同名 TXT。项目不生成 CSV。所有公开循环数均由候选独立求解到 `1e-6` 得到，选择阶段预测只用于筛选。
