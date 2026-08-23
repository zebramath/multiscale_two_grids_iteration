# two_grids_iteration 2.2.0

这是一个用于高对比扩散问题两网格粗空间机制研究的自包含 C++17
实验平台。项目重点不是提供生产级 AMG，而是在统一离散和两网格循环下，
可复现地比较局部化半径、局部求解精度、稀疏截断和非局部支撑选择。

## 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

若希望针对当前 CPU 编译，可增加
`-DTGI_NATIVE_OPTIMIZATION=ON`。所有实验默认把结果写入 `results/`。

## 实验程序

| 程序 | 内容 |
|---|---|
| `generate_test_fields` | 生成三类共享系数场 |
| `experiment1` | 固定几何支撑半径 |
| `experiment2` | 局部 PCG 容差 |
| `experiment3` | 全局高精度基的截断诊断 |
| `experiment4` | 原始固定预算强连接扩展 |
| `experiment5` | 固定预算、残差预算、强度归一化和扩展后压缩消融 |
| `experiment6` | 多右端项和误差传播谱代理 |
| `experiment7` | 种子、对比度和尺度扫描 |

新增实验 5–7 同时生成便于后处理的 CSV 文件。

## 新增支撑算法

### 固定强路径扩展

原始 `build_residual_strong_supports` 保留为稳定基线。实验 5 的优化版本只对
真正发生支撑变化的列进行温启动重求解，并使用更紧的修正容差，避免无变化列
重复求解，也避免温启动过早停止。

### 逐轮残差预算

`residual_budget_support.hpp` 提供探索性算法：

1. 计算当前插值的缩放细点残差；
2. 按未消除残差能量标记列；
3. 从残差种子沿强连接路径扩展一小批节点；
4. 只温启动重求解变化列；
5. 重新计算残差并按目标比例停止。

它支持全局最大值、逐行最大值和对称对角归一化三种强连接定义。当前数值结果
表明，全局残差预算可能把资源过度集中到少数列，因此它被保留为研究消融，
而不是默认优于固定预算的方法。

### 扩展后压缩

`support_pruning.hpp` 先在较大强连接支撑上求解，再按基函数幅值压缩非局部
节点，最后在保留支撑上重新能量极小化。该方法可在最终支撑接近 `K=64` 时，
获得比直接 `K=64` 更好的循环数；其额外 setup 成本是否值得取决于场拓扑和
右端项数量。

## 推荐运行方式

快速检查：

```bash
scripts/run_research_suite.sh quick
```

完整复现实验：

```bash
scripts/run_research_suite.sh full
```

单独做小型稳健性扫描：

```bash
./build/experiment7 --fine=64 --coarse=8 --seeds=3
```

目标尺度扫描示例：

```bash
./build/experiment7 --fine=128 --coarse=16 --seeds=5 --contrast=1e4
./build/experiment7 --fine=256 --coarse=16 --seeds=5 --contrast=1e6
```

实验 6 默认使用 6 个右端项和两个确定性谱估计起点。日常测试可降低成本：

```bash
./build/experiment6 --spectral-iters=30 --max-cycles=4000
```

## 结果解释边界

- `exact`/`ref` 是 PCG 高精度数值参考，不是符号意义的精确解。
- 两网格循环数依赖粗点集合、平滑器和右端项；实验 6 用于降低单一右端项偏差。
- 固定粗点下的全局能量基仍可能在通道问题上收敛很慢，支撑优化不能替代粗点
  选择或谱富集。
- 计时应在 Release 构建、固定线程、预热和多次重复条件下解释。

