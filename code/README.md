# Numerical implementation (v9.1)

代码只保留研究主线所需的组件：从几何插值出发的全局列 Jacobi--PCG 路径、固定路径点、energy endpoint、真实两网格谱诊断、有限传播诊断以及受控的三层扩展。旧版经验 adaptive API、oracle 扫描和 stopping-rule 消融已删除。

## 正式实验

| 程序 | 目的 |
|---|---|
| `experiment1_finite_path_comparison` | 13 个问题上分别比较 $m=(1/h)/4,(1/h)/3,(1/h)/2$ 与 endpoint |
| `experiment2_spectral_path` | 两种拓扑的 $m=1,\ldots,128$ 真实 $\rho_{TG}$ 与 RHS 因子全路径 |
| `experiment3_scaling_propagation` | fixed-$H$/fixed-$q$ 双尺度及有限图传播 |
| `experiment4_local_diagnostic_export` | 导出小规模精确矩阵与 PCG 路径；Python 完成 gap/方向诊断 |
| `experiment5_robustness` | 5 个 coefficient seed、6 个 RHS、5 次交替计时 |
| `experiment6_endpoint_comparison` | 几何插值与高精度 energy endpoint 的受控基线 |
| `experiment7_multilevel_pilot` | 两个三层 V-cycle 可行性案例 |

`scripts/analyze_local_diagnostic.py` 直接构造 $A^{1/2},B^{1/2},S^{-1/2},T_F,T_*,Z_m$，并交叉核对图空间公式与稠密特征值。`spectral_diagnostics.hpp` 对 $A$-自伴正半定的对称两网格误差传播算子做矩阵自由 Lanczos；CSV 中的 stage difference 是半程与全程 Ritz 值之差。

## 构建与运行

```bash
./scripts/run_all_experiments.sh quick
./scripts/run_all_experiments.sh full
```

脚本以 `-O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror` 构建；检测到 OpenMP 时启用。正式 full 模式依次生成全部文本、CSV 和 PNG。谱路径支持 `--topology=cross|ring|both`，并提供 `--endpoint-only` 以便长运行中断后只恢复端点，不必重跑有效路径。

## 数值口径

- PDE：单位方形二维五点扩散，Dirichlet 边界。
- 粗点：规则 C/F 划分；$q=H/h$。
- 插值：粗点单位注入；F 行从双线性几何初值沿列 PCG 更新。
- 磨光：一次前向和一次后向 Gauss--Seidel。
- 粗解：精确稀疏 Cholesky；$A_c=P^TAP$。
- 线性求解终止：相对残量 $10^{-6}$，通常最多 20,000 cycles。
- energy endpoint：每列相对残量不超过 $10^{-10}$（局部诊断为 $10^{-12}$）。
- 固定分数步数统一使用最近整数函数 `fixed_path_steps`。

所有正式结果位于 `results/`，其索引说明见 `results/README.md`。
