# multiscale_two_grids_iteration v8.4 验证记录

验证日期：2026-09-07  
代码版本：8.4.0  
正式运行线程数：4

## 文件完整性

用户提供的理论主文件与三份样式库按原始字节纳入项目，未作内容修改。

| 文件 | SHA-256 |
|---|---|
| `theory.tex` | `43b6477fa7d684973d0a1bbcabdaf564e422f7f41fd1f0023625fcd00a1cae1a` |
| `packages.sty` | `c10ff69f4d6e81b98491dfc7917ae2aae3148e7e50ea07273836f1839209670c` |
| `mathstyle.sty` | `1a9032419c3a6b37fb395f8533fcae0707a4aede0ff52ee72543e3c8e70abaf1` |
| `reportstyle.sty` | `ecaa9ad7a53709a5bcf8df1932df07298cd60984a04673b7efd287e34ebc75ad` |

`theory.tex` 与本次用户附件逐字节一致；三份样式库沿用既有附件，未作修改。理论文件在
替换后仅作读取与完整性核验，不再改动。

## 构建与静态检查

- 编译器：GCC 13.3.0，C++17，`-O3 -DNDEBUG`。
- 启用 `-Wall -Wextra -Wpedantic -Werror`、shadow、conversion、sign-conversion、
  duplicated-condition、logical、null-dereference 和 format 检查；八个实验入口全部通过。
- `run_validation.sh` 通过 shell 语法检查，`plot_path_scan.py` 通过 Python 字节码检查。
- quick 回归的三个代表问题全部收敛，循环数与既有正式结果逐项一致。
- AddressSanitizer 与 UndefinedBehaviorSanitizer 定向检查覆盖零右端项、零循环预算、非法
  容差、PCG 路径回退拒绝和常规两网格收敛，未报告地址或未定义行为错误。

## 数值结果完整性

Exp1--3 与 Exp5--7 保留 v8.1.0 的完整重算结果，Exp8 保留 v8.3.0 的新增结果。v8.4
未重复运行这些实验；Exp1--3、Exp5--8 的数据文件与 v8.3 逐字节一致。Exp4 由 v8.4
定向复算。

| 实验 | 完整性检查 |
|---|---|
| Exp1 | 13 个问题，adaptive/global-reference 均 13/13 收敛，循环和 4459/28972 |
| Exp2 | 两个拓扑各 128 点，`m=1,...,128` 连续无缺失；两条原始能量序列逐步严格下降 |
| Exp3 | 7 个设计问题与 3 个冻结验证问题全部完成 |
| Exp4 | 5 个 coefficient seed、6 个 RHS、双方预热及交替顺序的 5 次正式计时全部完成 |
| Exp5 | 6 个问题、3 种停止策略全部完成 |
| Exp6 | 3 层固定物理加密完成，共享节点系数失配数为 0 |
| Exp7 | 2 个三层层次、2 种插值方法全部完成 |

## Exp4 定向复算

Exp4 由 v8.4.0 重新运行。五个 coefficient seed 与六个 RHS 的循环数、收敛状态和有效
因子均与原结果一致；计时表按既定呈现要求只保留五次正式测量的算术平均值。

| 方法 | setup/ms | solve/ms | total/ms | 循环数 |
|---|---:|---:|---:|---:|
| adaptive | 391.970 | 229.785 | 621.755 | 242 |
| global-reference | 2264.647 | 7489.791 | 9754.438 | 3227 |

## 中心问题插值对照

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

## 理论稿编译检查

已调用 XeLaTeX（TeX Live 2023）实际编译 `theory.tex`。当前环境缺少 `ctexart.cls`，
编译在文档类加载阶段停止，尚未读取正文或自定义样式；因此无法在本环境生成 PDF。
只读结构检查确认 159 组 LaTeX 环境全部配对、160 个标签无重复、107 个 `\cref` 目标
均存在。理论 TEX 和样式库保持用户提供的原始字节，不为规避环境依赖而修改。安装完整
TeX Live 中文组件后可在项目根目录运行：

```bash
xelatex theory.tex
xelatex theory.tex
```
