# two_grids_iteration v3.0.0

本项目研究二维高对比扩散问题中，能量插值的空间局部化、有限步 PCG
局部化和稀疏支撑扩展。v3.0.0 在 v2.11.0 的六组基线实验上增加：

- 可在检查点间继续推进、不会重复 Krylov 计算的全局 PCG 路径；
- 同时使用两网格谱探针和代表性右端项尾窗收缩率的自适应停止器；
- 带全局预算、公平性下限和实测停滞判据的 `frontier-gain` 支撑扩展；
- 交叉、弯曲、对角和平行通道四种拓扑；
- 多随机种子、多对比度和多组 `H/h` 的 18 问题鲁棒性矩阵；
- 增量路径等价性、自适应决策和支撑预算回归测试。

## 构建

需要 C++17 编译器、CMake 3.16 以上版本和线程库；OpenMP 可选。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

如果环境没有 CMake，单个程序也可直接编译：

```bash
g++ -std=c++17 -O2 -pthread -fopenmp -Isrc \
    experiments/experiment9.cpp -o experiment9
```

## 验证

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

`quick` 运行单元测试、两个核心消融实验和 6 问题代表性扫描；`full`
运行 128×128 原始通道基准、四拓扑支撑对比和完整 18 问题扫描。
线程数可通过 `TGI_THREADS` 设置。

## v3 实验

- `experiment7`：固定 PCG 检查点与自适应选择轨迹；
- `experiment8`：v2 预算式支撑与 v3 前沿收益支撑的消融；
- `experiment9`：网格、`H/h`、对比度、种子和通道拓扑鲁棒性矩阵。

结果写入 `results/*.txt` 和 `results/*.csv`。计时是单次墙钟时间，适合在同一
机器内横向比较，不应跨机器直接比较绝对数值。

## 主要接口

- `multigrid/global_pcg_path.hpp`：`GlobalEnergyPcgPath`；
- `multigrid/adaptive_global_pcg.hpp`：
  `build_adaptive_global_pcg_interpolation`；
- `multigrid/frontier_gain_support.hpp`：
  `build_frontier_gain_interpolation`。

自适应 PCG 可传入代表性右端项。若求解任务不存在代表性右端项，应保留多个谱
探针，并提高 `power_iterations`；不能用普通 PCG 方程残差替代两网格质量探针。

## 结论边界

代码验证范围是二维结构网格、固定几何粗点、一次前后向 Gauss--Seidel 和
Galerkin 粗校正。自适应停止是有限代价的模型选择规则，不是任意 SPD 系统上的
全局最优保证；`frontier-gain` 的优势主要表现为相同迭代数下更低的插值密度，
其多轮构造时间可能高于旧策略。
