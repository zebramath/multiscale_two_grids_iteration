# Numerical implementation

代码实现从双线性几何插值出发的全局列 Jacobi--PCG 路径、固定路径点、energy endpoint、真实两网格谱诊断、有限传播诊断和三层 V-cycle。

## 实验

| 程序 | 目的 |
|---|---|
| `experiment1_finite_path_comparison` | 13 个问题上比较 $m=(1/h)/4,(1/h)/3,(1/h)/2$ 与 endpoint |
| `experiment2_spectral_path` | 两种拓扑的 $m=1,\ldots,128$ 真实 $\rho_{TG}$ 与 RHS 因子路径 |
| `experiment3_scaling_propagation` | fixed-$H$/fixed-$q$ 双尺度及有限图传播 |
| `experiment4_local_diagnostic_export` | 导出小规模矩阵与 PCG 路径，供 gap/方向诊断使用 |
| `experiment5_robustness` | 5 个 coefficient seed、6 个 RHS、5 次交替计时 |
| `experiment6_endpoint_comparison` | 几何插值与 energy endpoint 基线 |
| `experiment7_multilevel_pilot` | 两个三层 V-cycle 案例 |

`scripts/analyze_local_diagnostic.py` 构造 $A^{1/2},B^{1/2},S^{-1/2},T_F,T_*,Z_m$，并核对图空间公式与稠密特征值。`spectral_diagnostics.hpp` 对 $A$-自伴正半定的对称两网格误差传播算子执行矩阵自由 Lanczos；CSV 中的 stage difference 为半程与全程 Ritz 值之差。

## 构建与运行

```bash
./scripts/run_all_experiments.sh quick
./scripts/run_all_experiments.sh full
```

脚本使用 `-O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror`。OpenMP 可用时启用。谱路径接受 `--topology=cross|ring|both`。

## 数值口径

- PDE：单位方形二维五点扩散，Dirichlet 边界。
- 粗点：规则 C/F 划分，$q=H/h$。
- 插值：粗点单位注入；F 行从双线性几何初值沿列 PCG 更新。
- 磨光：一次前向和一次后向 Gauss--Seidel。
- 粗解：精确稀疏 Cholesky，$A_c=P^TAP$。
- 线性求解相对残量：$10^{-6}$。
- energy endpoint 列相对残量：$10^{-10}$；局部诊断为 $10^{-12}$。
- 固定分数步数使用最近整数函数 `fixed_path_steps`。

正式结果位于 `results/`。
