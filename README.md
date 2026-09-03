# multiscale_two_grids_iteration v7.7

本项目研究高对比扩散问题中的有限 Krylov 能量插值。固定粗点与磨光器后，能量极小化
路径和两网格收敛目标可以呈现不同的变化规律：插值能量沿 PCG 路径持续下降，而两网格
收敛因子可能在有限步达到更小值。项目据此分析有限步插值，并构造自适应有限步规则。

## 方法

记每个坐标方向的细网格区间数为 $n=1/h$，粗网格区间数为 $n_H=1/H$。`adaptive`
根据 $n_H$ 与矩阵对角尺度比选择

$$
m\in\{\operatorname{round}(n/8),\operatorname{round}(n/4),
       \operatorname{round}(n/3),\operatorname{round}(n/2)\}.
$$

标准 PCG 条件数估计给出达到固定能量衰减所需的 $m=O(h^{-1})$ 充分尺度；候选比例与
尺度比阈值由设计问题组标定。在线阶段确定单个候选，执行一次 PCG 路径推进、插值组装
和 Galerkin 粗算子构造。

## 主要内容

- PCG 正交能量账本与归一化插值误差；
- Galerkin 粗矩阵、粗空间主角和两网格投影的统一表达；
- 重主导奇异值下的一阶有利方向、局部 Pareto 律和有限步反单调条件；
- 能量误差到两网格性能的稳定界与递归 V-cycle 误差分解；
- 四个数值主题、七个独立实验，覆盖两网格主结果与逐步扫描、自适应规则评估、固定物理
  场加密和三层 V-cycle pilot。

当前分析采用二维规则网格、固定粗点、对称 Gauss--Seidel 和 Galerkin 粗算子。固定物理
场试验和三层试验分别考察标准网格加密与递归构造；更一般的离散和层次结构是后续方向。

## 文件

| 文件或目录 | 内容 |
|---|---|
| `research_report.md` | 研究定位、实验设计、完整数值结果和结论 |
| `theory.tex` | 定理、命题、证明、假设和参考文献 |
| `code/` | C++17 实现、七个实验入口和运行脚本 |
| `code/results/` | v7.7 完整复跑结果、逐步扫描数据和曲线 |
| `VALIDATION_v7.7.md` | 源码精简、构建、文稿和完整运行验证记录 |

快速验证：

```bash
cd code
./scripts/run_validation.sh quick
```

完整复现：

```bash
cd code
./scripts/run_validation.sh full
```

理论稿使用 `ctexart`，建议以 XeLaTeX 和完整的 TeX Live 中文组件编译。
