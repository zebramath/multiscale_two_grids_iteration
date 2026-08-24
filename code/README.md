# two_grids_iteration v3.7.0

## 目录结构

`multigrid` 只保留四个核心模块：

| 模块 | 内容 |
|---|---|
| `energy_interpolation.hpp` | 几何、局部、全局能量插值，强连接距离插值和参考剪枝 |
| `support_expansion.hpp` | 强连接支撑扩张和残差预算自适应支撑 |
| `global_pcg.hpp` | 增量全局 PCG 路径和自适应检查点选择 |
| `two_grid_solver.hpp` | Galerkin 粗算子、对称平滑、粗解和两网格循环 |

v3.7 的 `CoarseSetupReport` 还直接记录
$\operatorname{tr}(P^TAP)$、$C_P$、两层 operator complexity、对称化 Cholesky
填充比。experiment2 使用这些量和 60 步两网格幂迭代谱代理解释有限 PCG 路径；
这些诊断不计入公开 Setup/Total 计时。

`experiments` 只包含六个可执行入口。原有测试角度位于 `src/experiment/studies`，由对应入口统一组织并写入一份 TXT，不再分别形成零散实验。

## 自适应算法

选择器复用一条增量式全局能量 PCG 路径，最多建立四个层次：几何初值、最小步锚点、区间中点和一个补充点。每个候选运行至多 24 个两网格 pilot，不在 setup 阶段求解到完整容差。候选位置不读取系数场名称、对比度、种子或预设最优步数。

几何初值或最小步锚点已经足够好时立即停止；否则检查中点。中点仍较难时，仅用归一化 PCG 能量残差与求解容差的平方根比较，选择锚点—中点内点或中点—上界前向点。在最小预测的 13% 范围内选择最小 `m`。选中的插值矩阵、Galerkin 粗算子、排序和稀疏 Cholesky 因子直接交给正式求解与多右端项工作负载复用。

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

## 六个整合实验

| 入口 | 内部测试角度 |
|---|---|
| `experiment1_localization` | 支撑半径、局部 PCG 容差、全局剪枝、支撑扩张 |
| `experiment2_pcg_path` | 有限 PCG 预算、步长 2 密集扫描、能量/复杂度/谱诊断、pilot 排序 |
| `experiment3_adaptive_oracle` | 主问题自适应轨迹、八问题 step-2 oracle、pilot×松弛初步消融 |
| `experiment4_robustness` | 18 问题网格、对比度、种子与六类通道拓扑矩阵 |
| `experiment5_diagnostics` | 八问题附加诊断和八问题压力测试 |
| `experiment6_workload` | 五右端项摊销与 break-even |

每个入口只生成同名 TXT。项目不生成 CSV。所有公开循环数均由选中候选独立求解到 `1e-6` 得到，选择阶段预测只用于筛选。
