# two_grids_iteration v3.2.0

本项目研究二维高对比扩散问题中，能量插值、有限步 PCG 与两网格收敛之间的
非单调关系。v3.2.0 将自适应器的目标从“尽量命中最小循环数”改为“在循环数
损失可接受时显著压缩选择 Setup”。

## v3.2 的低预算选择器

默认策略只做以下工作：

1. 对几何插值运行短 pilot；预测循环数已经很小时直接接受；
2. 沿同一个 `GlobalEnergyPcgPath` 依次产生少量检查点；
3. 一般使用 `m_min`、`m_min+12` 和 `m=40`，大网格省略中间点；
4. `H/h` 较大时适当延长 pilot，减轻慢瞬态造成的预测偏差；
5. 不再把两个最终候选都运行到真实容差，真实循环数在独立评价阶段测量。

`AdaptiveGlobalPcgOptions::cost_aware_mode=false` 可切回 v3.1 的
“粗筛—回溯—加密—真实确认”高质量模式。v3.2 修复了该模式中的重复构造：
候选两网格对象在 pilot 和 confirmation 之间复用。

## 主要实测结论

- 原始 `128 x 128` 交叉通道：v3.2 选择 `m=40`，234次循环；离步长2
  oracle 的231次仅差1.3%。同一进程中低预算 Setup 约496 ms，高质量 staged
  模式约1495 ms，下降约66.8%。
- 18问题矩阵：低预算模式总循环数1881，高质量模式1813，只增加3.75%；
  低预算模式18/18均不超过高质量模式的1.3倍。
- 同一18问题矩阵：低预算 Setup 总和约1052 ms，高质量模式约1926 ms，
  下降约45.4%，18/18逐例更低。
- 六问题步长2 oracle：平均差距约6.85%，最大20%，3/6精确命中。
- 最终冻结验证的8个新种子/新系数场中，7例精确命中步长4 oracle；唯一困难
  交叉通道差距33.8%。

墙钟时间是单机单次或少量重复测量，只适合同机比较。完整数字见 `results/` 和
研究方案 v3.2。

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

线程数通过 `TGI_THREADS` 设置。脚本会在每个构建、测试和实验步骤开始、结束时输出
状态；experiment9、11--14 还会逐案例输出进度。每个步骤默认最多运行900秒，必要时
可用 `TGI_STEP_TIMEOUT_SECONDS` 调整，例如：

```bash
TGI_STEP_TIMEOUT_SECONDS=1800 ./scripts/run_validation.sh full
```

核心两网格循环、局部 PCG、自适应候选和确认过程均有显式整数上界。源码中的
`while (true)` 仅用于有限原子任务队列：索引达到预先确定的列数或支撑数后退出。

## v3.2 实验

- `experiment7`：原始通道、固定候选和低预算轨迹；
- `experiment9`：18问题矩阵以及低预算/高质量模式消融；
- `experiment10`：不同 pilot 长度的残差诊断；
- `experiment11`：六问题步长2 oracle；
- `experiment12`：额外种子暴露失败后得到的诊断集；
- `experiment13`：策略冻结后的新种子和非通道系数场验证；
- `experiment14`：五种右端项、Setup 摊销和盈亏平衡实验。

## 主要接口

- `multigrid/global_pcg_path.hpp`：可继续推进并保存检查点的 PCG 路径；
- `multigrid/adaptive_global_pcg.hpp`：低预算和高质量两种选择模式；
- `build_adaptive_global_pcg_interpolation(...)`：统一入口。

低预算模式报告中的 `estimated_selected_cycles` 是 pilot 预测，不是确认值；
`selected_cycles_confirmed=false` 明确标识这一点。实验表中的 `Cycles` 均由选中
候选重新求解到真实容差得到，不把预测值冒充实测值。

## 结论边界

当前范围仍是二维结构网格、固定几何粗点、一次前后向 Gauss--Seidel 与 Galerkin
粗校正。默认候选位置和 pilot 长度是实验驱动参数。单个代表右端项上的选择不保证
最坏情形谱半径最优，也不保证任意 SPD 系统上都在 oracle 的30%以内。

若事先已经知道固定 `m=40` 合适，而自适应器最终也选择 `m=40`，自适应器无法
回收额外选择成本；它的价值来自候选事先未知、几何方法可能失败或需要复用多个
右端项的场景。
