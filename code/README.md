# two_grids_iteration v3.1.0

本项目研究二维高对比扩散问题中，能量插值、有限步 PCG 与两网格收敛之间的
非单调关系。v3.1.0 在 v3.0.0 的增量 PCG 路径上重点重构自适应停止器：

- 用“已完成循环 + 预计剩余循环”代替只看尾部收缩率；
- 采用粗筛、向前期回溯、步长2加密和少量真实容差确认；
- 缓存粗筛路径上的细检查点，不再为加密阶段重跑全局 PCG；
- 对实际低循环候选提供确认后快速接受；
- 对中等候选只增加一个前瞻候选，限制选择开销；
- 保留交叉、弯曲、对角和平行通道及18问题鲁棒性矩阵；
- 增加代表性问题的离线步长2 oracle 对比。

实验表明收益有限的 `frontier-gain` 支撑扩展已经从 v3.1.0 删除。原来的
`adaptive-budget-v2` 仍保留在 `residual_budget_support.hpp` 和 experiment4 中，
但不是 v3.1 的研究主线。

## 构建

需要 C++17 编译器、CMake 3.16 以上版本和线程库；OpenMP 可选。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

没有 CMake 时可直接编译：

```bash
g++ -std=c++17 -O2 -pthread -fopenmp -Isrc \
    experiments/experiment9.cpp -o experiment9
```

## 验证

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

`quick` 运行测试、缩小版鲁棒性矩阵和 oracle 对比；`full` 还会运行
128×128 原始通道及完整18问题矩阵。线程数通过 `TGI_THREADS` 设置。

## v3.1 实验

- `experiment7`：原始通道上的固定候选与分阶段自适应轨迹；
- `experiment9`：多网格、`H/h`、对比度、种子和拓扑鲁棒性矩阵；
- `experiment10`：PCG 步数与不同试运行长度的残差诊断；
- `experiment11`：六个代表问题与离线步长2 oracle 的差距。

结果写入 `results/*.txt` 和 `results/*.csv`。墙钟时间是单次同机测量，跨机器
不可直接比较。

## 主要接口

- `multigrid/global_pcg_path.hpp`：`GlobalEnergyPcgPath`；
- `multigrid/adaptive_global_pcg.hpp`：
  `build_adaptive_global_pcg_interpolation`。

v3.1 的最终确认针对传入的代表性右端项，因此调用者必须提供该右端项。算法优化
的是这一工作负载下的循环数；对一组性质不同的右端项，应扩展为多代表向量确认，
不能直接宣称得到最坏情形谱半径最优。

## 结论边界

当前范围仍是二维结构网格、固定几何粗点、一次前后向 Gauss--Seidel 与 Galerkin
粗校正。分阶段搜索能显著接近离线候选最优并降低平均选择成本，但不保证任意 SPD
问题上的全局最优 PCG 步数。
