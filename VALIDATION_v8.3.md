# multiscale_two_grids_iteration v8.3 验证记录

验证日期：2026-09-06  
代码版本：8.3.0  
新增实验线程数：4

## 文件完整性

用户提供的理论主文件与三份样式库按原始字节纳入项目，未作内容修改。

| 文件 | SHA-256 |
|---|---|
| `theory.tex` | `a0597dcb49e50e4a71f15fb4420fcfd2b82153f120186998d81b2627a667fd8b` |
| `packages.sty` | `c10ff69f4d6e81b98491dfc7917ae2aae3148e7e50ea07273836f1839209670c` |
| `mathstyle.sty` | `1a9032419c3a6b37fb395f8533fcae0707a4aede0ff52ee72543e3c8e70abaf1` |
| `reportstyle.sty` | `ecaa9ad7a53709a5bcf8df1932df07298cd60984a04673b7efd287e34ebc75ad` |

这些哈希与四个用户附件的原始哈希逐项一致。

## 构建与静态检查

- 编译器：GCC 13.3.0，C++17，`-O3 -DNDEBUG`。
- 启用 `-Wall -Wextra -Wpedantic -Werror`、shadow、conversion、sign-conversion、
  duplicated-condition、logical、null-dereference 和 format 检查；八个实验入口全部通过。
- `run_validation.sh` 通过 `sh -n`。
- `plot_path_scan.py` 已在 v8.1 完整复算中通过 Python 3.12.13 字节码检查，绘图使用
  Matplotlib 3.10.8；v8.3 未改动该脚本。
- v8.3 仅执行新增 Exp8；其余实验结果文件经逐文件比较保持不变。

## 原有数值结果

Exp1--7 保留 v8.1.0 可执行文件的完整重算结果，v8.3 未重复运行这些实验。

| 实验 | 完整性检查 |
|---|---|
| Exp1 | 13 个问题，adaptive/global-reference 均 13/13 收敛，循环和 4459/28972 |
| Exp2 | 两个拓扑各 128 点，`m=1,...,128` 连续无缺失；两条原始能量序列逐步严格下降 |
| Exp3 | 7 个设计问题与 3 个冻结验证问题全部完成 |
| Exp4 | 5 个 coefficient seed、6 个 RHS、双方预热及交替顺序的 5 次正式计时全部完成 |
| Exp5 | 6 个问题、3 种停止策略全部完成 |
| Exp6 | 3 层固定物理加密完成，共享节点系数失配数为 0 |
| Exp7 | 2 个三层层次、2 种插值方法全部完成 |

## v8.3 新增中心问题对照

Exp8 采用 128/16 网格、对比度 $10^4$、cross-channel、seed 1 和常数右端项。几何插值
与能量极小插值共用矩阵、粗点、Galerkin 构造、对称 Gauss--Seidel、零初值及
$10^{-6}$ 相对残量停止标准，仅插值矩阵不同。两种方法均单次运行至收敛，不设循环数
上限或脚本超时。

| 插值 | 循环数 | 最终相对残量 | 有效收敛因子 |
|---|---:|---:|---:|
| 几何插值 | 85524 | $9.99896215\times10^{-7}$ | 0.999838472 |
| 能量极小插值 | 3227 | $9.97423567\times10^{-7}$ | 0.995727131 |

几何插值循环数是能量极小插值的 26.503 倍。新增结果保存于
`code/results/experiment8_interpolation_endpoint.txt`。

Exp2 自动核验结果：

| 拓扑 | 扫描区间内最小值 | adaptive 点 | 最大相邻能量增量 | 末点归一化能量差 |
|---|---:|---:|---:|---:|
| cross-channel | `m=38`, 0.941763675 | `m=43`, 0.944406934 | -29.9540004730 | 1.233657970477e-08 |
| winding-ring | `m=53`, 0.941012129 | `m=43`, 0.943169306 | -0.0119991302 | 5.810126337407e-12 |

CSV 舍入前的最大相邻能量增量记录于 `experiment2_step_scan.txt`；两者均为负。两张
双面板 PNG 已人工检查，能量与 $\rho_{\mathrm{eff}}$ 面板、最小值和 adaptive 标记均正常。

## 定向内存安全复验

首次运行新的 warm-start timing 分支时，AddressSanitizer 定位到一个临时 prolongation
对象的悬空引用。代码已改为保存具名 `SparseMatrix` 后再构造 `TwoGridCycle`。修复后 Exp4
在 AddressSanitizer 下完整运行并以 0 退出，所有 seed、RHS、预热和 5 次计时均通过，未
报告地址、越界或生命周期错误。由于 LeakSanitizer 在受监控环境中不可用，终态复验仅
关闭 leak 检测，保留全部 AddressSanitizer 地址检查。

## 理论稿编译检查

已调用 XeLaTeX（TeX Live 2023）实际编译 `theory.tex`。当前环境缺少 `ctexart.cls`，
编译在文档类加载阶段停止，尚未读取正文或自定义样式；因此无法在本环境生成 PDF。
理论 TEX 和样式库保持用户提供的原始字节，不为规避环境依赖而修改。安装完整 TeX Live
中文组件后可在项目根目录运行：

```bash
xelatex theory.tex
xelatex theory.tex
```
