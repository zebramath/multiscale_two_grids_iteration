# two_grids_iteration v3.3.0

## 当前算法

v3.3 使用同一条增量式全局能量 PCG 路径构造候选插值。选择器最多建立四个两网格层次：

1. 几何插值 `m=0`；
2. 最小正检查点 `m_min`；
3. 由前两个短 pilot 的对数收敛斜率投影出的检查点，投影上界为剩余区间中点；
4. 仅在预测仍很差或 PCG 能量残差显示已越过快速衰减区时，向前探测或回到局部中点。

候选位置不依赖系数场名称、对比度、种子或预设的 `m=40`。每个候选只运行24个两网格 pilot，不运行到真实求解容差。尾部模型使用重叠窗口中最慢的对数衰减率，并记录中位数绝对偏差给出的不确定度。选择时在最佳预测的8%松弛范围内取最小 `m`。

PCG 路径同时报告最大相对残差、RMS 相对残差和

\[
\eta_E=\left(\frac{\sum_j r_{j,m}^{T}D^{-1}r_{j,m}}
{\sum_j r_{j,0}^{T}D^{-1}r_{j,0}}\right)^{1/2},
\]

后者与理论中的残差对偶能量指标直接对应，用于决定第四个候选的方向。

## 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

没有 CMake 时可直接编译：

```bash
g++ -std=c++17 -O3 -DNDEBUG -pthread -fopenmp -Wall -Wextra -Wpedantic \
    -Isrc experiments/experiment7.cpp -o experiment7
```

## 验证

```bash
./scripts/run_validation.sh quick
./scripts/run_validation.sh full
```

可通过 `TGI_THREADS` 和 `TGI_STEP_TIMEOUT_SECONDS` 调整线程数与单步超时。

## 实验索引

| 脚本 | 研究角度 |
|---|---|
| experiment1 | 空间支撑半径与全局能量参考 |
| experiment2 | 固定支撑上的局部 PCG 容差 |
| experiment3 | 全局插值的幅值剪枝扰动 |
| experiment4 | 固定与残差预算支撑扩张 |
| experiment5 | 全局 Jacobi、有限 PCG 和精确能量插值 |
| experiment6 | 通道问题的步长2密集扫描与独立粗解交叉验证 |
| experiment7 | 原始大通道问题上的 v3.3 选择轨迹 |
| experiment9 | 18问题鲁棒性矩阵 |
| experiment10 | pilot 残差排序诊断 |
| experiment11 | 六问题步长2 oracle |
| experiment12 | 额外种子和系数场诊断 |
| experiment13 | 额外策略压力测试 |
| experiment14 | 五右端项摊销与 break-even |

所有公开循环数均由选中候选独立求解到 `1e-6` 得到；`estimated_selected_cycles` 仅是选择阶段预测。输入尺寸与选项检查已从计算路径移除，调用者应保证输入合法；SPD 分解失败、迭代不收敛等数学有效性检查仍保留。
