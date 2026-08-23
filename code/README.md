# two_grids_iteration v3.4.0

## 当前算法

项目使用同一条增量式全局能量 PCG 路径构造候选插值。自适应器最多建立四个两网格层次：

1. 几何插值 `m=0`；
2. 最小正检查点 `m_min`；
3. `[m_min,m_max]` 的量化中点；
4. 根据中点 pilot 预测和归一化 PCG 能量残差，至多向前或向后二分一次。

候选位置不读取系数场名称、对比度、种子或预设最优步数。每个候选只运行 24 个两网格 pilot，不在 setup 阶段求解到真实容差。若中点已经足够好则立即停止；只有中点仍很差时才添加一次细分。选择器在最佳预测的 8% 范围内取最小 `m`。

PCG 路径报告归一化预条件能量残差

\[
\eta_E=\left(\frac{\sum_j r_{j,m}^{T}D^{-1}r_{j,m}}
{\sum_j r_{j,0}^{T}D^{-1}r_{j,0}}\right)^{1/2}.
\]

向后细分阈值使用 `sqrt(solve_tolerance)`，不引入系数场相关参数。若所有 pilot 预测都达到硬上限，`eta_E` 作为唯一裁决量。

## 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
```

没有 CMake 时可直接编译：

```bash
g++ -std=c++17 -O3 -DNDEBUG -pthread -fopenmp \
    -Wall -Wextra -Wpedantic -Werror -Isrc \
    experiments/experiment7.cpp -o experiment7
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
| experiment5 | 几何初值、有限全局 PCG 和精确全局能量插值 |
| experiment6 | 通道问题的步长 2 密集扫描与独立粗解交叉验证 |
| experiment7 | 原始大通道问题上的自适应选择轨迹 |
| experiment9 | 18 问题鲁棒性矩阵 |
| experiment10 | pilot 残差排序诊断 |
| experiment11 | 六问题步长 2 oracle |
| experiment12 | 额外种子和系数场诊断 |
| experiment13 | 策略压力测试 |
| experiment14 | 五右端项摊销与 break-even |

所有公开循环数均由选中候选独立求解到 `1e-6` 得到；选择阶段预测只用于筛选。调用者负责输入尺寸和选项的合法性；SPD 分解失败、PCG breakdown 和显式要求收敛但未收敛等数学有效性检查仍保留。
