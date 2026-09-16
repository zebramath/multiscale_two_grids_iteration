# Finite PCG interpolation for high-contrast two-grid solvers — v10.0

本项目研究高对比、多尺度扩散离散所产生的困难 SPD 线性方程组。固定粗点、单位注入和对称 Gauss--Seidel 磨光后，从几何插值出发，用逐列 Jacobi--PCG 逼近 ideal/energy-minimizing interpolation。核心发现是：PCG 每个非零步都严格降低插值能量，但真实两网格谱半径并不随之单调降低；有限步中间点可以显著优于高精度能量极小端点。

## 核心结果

- 归一化图坐标 \(Z=A_{FF}^{1/2}(W-W^*)S^{-1/2}\) 同时给出精确能量差、Galerkin 粗矩阵和到 ideal 粗空间的主角。
- 对一次前向、一次后向 Gauss--Seidel 与精确粗解，真实对称两网格谱半径满足
  \[
  \rho_{TG}(Z)=\left\|(I+ZZ^T)^{-1/2}(T_F-ZT_*)\right\|_2^2.
  \]
- 在端点最大奇异值简单的条件下，论文给出带显式二阶余项的局部展开、PCG 相邻步变好/变坏的双向判据，以及任意路径位置的谱认证。
- PCG 校正满足严格的有限图传播：第 \(m\) 步校正支撑不超过初始残量的半径 \(m-1\) 邻域。
- 二维规则粗化下，Jacobi 条件数由统一容量尺度 \(\Lambda_{h,H}=h^2+[q^2\log(2q)]^{-1}\) 控制；固定物理粗网格 \(H\) 时自然 PCG 尺度为 \(\Theta(h^{-1})\)，固定粗化比 \(q=H/h\) 时为 \(O(1)\)（固定对比度）。

中心问题 \((1/h,1/H,\chi)=(128,16,10^4)\) 的循环数为：

| 插值 | 循环数 | 有效因子 | 插值密度 |
|---|---:|---:|---:|
| 几何双线性 | 85,524 | 0.999838472 | 1.3950% |
| 高精度能量极小端点 | 3,227 | 0.995727131 | 98.2099% |
| 有限 PCG，\(m=43\) | 242 | 0.944406934 | 29.4678% |

因此能量端点已经相对几何插值改善 26.50 倍，而有限 PCG 又在这一强基线上改善 13.33 倍。五次交替顺序计时的平均总时间为 0.679 s 对 10.193 s，即 15.01 倍。

## v10.0 相对 v9.4

- 审核了完整 72 次 Git 提交的文件树、删除记录和主要差异。
- 恢复了真正增强论文结论的两项历史内容：六类右端项稳健性，以及三层 V-cycle 可行性试验。
- 未恢复已被 v9 理论替代的旧自动选步/oracle、早期支撑剪枝路线、重复验证报告和过期结果。
- 将三个 README 合并为本文件；删除旧研究报告和编辑器配置。
- 将固定迭代求解循环抽象为可供两网格和多层共同使用的模板，避免复制求解状态逻辑。
- 所有 C++、Python、Shell 和 CMake 源文件已删除空白行，并以严格编译警告重新构建。
- 使用 Journal of Scientific Computing 官方 `svjour3` 模板和 `smallextended` 选项形成完整英文论文。
- 全量重跑 7 个实验入口并重新生成所有正式结果和图。

早期历史中有价值但不应重新并入主线的内容包括：旧 adaptive/oracle 选择器、停止阈值消融、支撑扩张/剪枝方法、旧版固定物理场表格和多个版本验证文档。其中固定物理场与 fixed-\(H\)/fixed-\(q\) 问题已由当前 experiment 3 更系统地覆盖；旧自动规则则与 v9.4 起采用的固定有限检查点定位不一致。

## 目录

| 路径 | 内容 |
|---|---|
| `paper/main.tex` | JSC 官方格式英文论文源文件 |
| `paper/main.pdf` | 编译后的论文 |
| `paper/references.bib` | 论文引用数据库 |
| `paper/svjour3.cls`, `paper/svglov3.clo`, `paper/spmpsci.bst` | JSC 官方下载包中的类与参考文献样式 |
| `theory.tex` | v9.4 中文详细推导底稿；论文已提炼其主定理和证明 |
| `code/src/` | 稀疏线性代数、扩散离散、PCG、两网格、多层与谱诊断 |
| `code/experiments/` | 7 个可复现实验入口 |
| `code/results/` | 全量重跑生成的正式文本、CSV 和图 |
| `code/scripts/run_all_experiments.sh` | 统一构建、运行、分析与绘图脚本 |

## 复现

依赖：支持 C++17 的编译器、POSIX shell、Python 3、NumPy、Matplotlib；CMake 可选。脚本在缺少 CMake 时自动采用直接编译。

```bash
cd code
./scripts/run_all_experiments.sh quick
./scripts/run_all_experiments.sh full
```

完整运行会覆盖 `code/results/` 中的正式结果。默认使用 4 线程；可通过 `TGI_THREADS` 修改：

```bash
TGI_THREADS=8 ./scripts/run_all_experiments.sh full
```

论文编译：

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -c
```

JSC 当前官方投稿指南要求 Springer LaTeX macro package 的 `smallextended` 选项，并要求同时提交源文件、样式、图片与编译 PDF。本目录已经按这一要求组织。投稿前只需在 `paper/main.tex` 填写作者单位、电子邮件、资助与最终数据仓库链接。

## 数值范围与结论边界

已证明结论采用一般 SPD 分块系统或二维规则五点扩散的明确假设；数值主结论采用固定粗点、单位注入、一次前向/后向 Gauss--Seidel、Galerkin 粗算子和精确粗解。三层试验是可行性证据，不被表述为多层收敛定理。`rho_TG` 是真实对称两网格谱半径，`rho_eff` 是给定右端项的有限残量统计，两者在代码和论文中始终分开。
