# two_grids_iteration v2.7.0

本项目研究高对比多尺度扩散问题中，全局能量极小插值的空间、代数和稀疏三类局部化，以及它们对固定两网格求解器的影响。

## 目录

```text
src/core/          CSR 矩阵、向量运算和稀疏 Cholesky
src/pde/           规则网格、系数场和扩散算子装配
src/multigrid/     插值构造、支撑扩展、剪枝和两网格循环
src/experiment/    配置、问题工厂、候选、评估和报告模块
experiments/       experiment1.cpp 至 experiment5.cpp
tests/             核心单元测试
scripts/           一键构建与运行脚本
results/           默认 128/16、对比度 1e4 的完整 CSV/TXT 结果
```

## 五组实验

| 编号 | 文件 | 研究因素 |
|---|---|---|
| experiment1 | `experiments/experiment1.cpp` | 几何插值、2/3/4 层局部能量基和全局能量基 |
| experiment2 | `experiments/experiment2.cpp` | 固定四层支撑上的 PCG 容差 |
| experiment3 | `experiments/experiment3.cpp` | 全局基的列相对幅值剪枝 |
| experiment4 | `experiments/experiment4.cpp` | 固定层、强连接扩展和预算式自适应扩展 |
| experiment5 | `experiments/experiment5.cpp` | 以全局 F 系统为目标的有限步 Jacobi/PCG |

数值输出列只保留结构、成本和收敛指标，不输出残差列。`Rho avg` 是实际求解全过程的平均收敛因子。

## 构建和运行

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --target experiments unit_core -j
ctest --test-dir build --output-on-failure

./build/experiment1
./build/experiment2
./build/experiment3
./build/experiment4
./build/experiment5
```

也可以执行：

```bash
./scripts/run_all.sh
```

所有实验接受：

```text
--fine=N --coarse=N --threads=N --contrast=X --max-cycles=N
```

默认值为 `--fine=128 --coarse=16 --threads=4 --contrast=1e4 --max-cycles=40000`。

## 结果口径

- `P nnz`、`P density %`：插值矩阵非零元数与相对稠密矩阵的密度。
- `Ac nnz`：Galerkin 粗矩阵非零元数。
- `L nnz`：粗矩阵 Cholesky 因子非零元数。
- `Build ms`：插值构造时间。
- `TG setup ms`：Galerkin 乘积、排序和粗矩阵分解时间。
- `Setup ms = Build ms + TG setup ms`。
- `Total ms = Setup ms + Solve ms`（只含一次端到端求解）。
- `Cycles`：达到固定求解容差所需两网格循环数。
- `Rho avg`：初始和终止范数比在全部循环上的几何平均因子。
- `Build iters`：逐列局部求解的平均迭代数；固定步方法就是步数。
- `Converged`：是否在循环上限内达到求解容差。

完整理论、实验设计、全部数值表和结论见同批交付的 `research_report_v2.7.md`。
