# multiscale_two_grids_iteration v7.6 验证记录

数值源码版本：7.6.0

理论与研究报告版本：v7.6

## 1. 结构与文本审校

- 七个实验程序和七份逐行结果保持独立，研究报告按四个数值主题组织：Exp1--2、
  Exp3--5、Exp6 和 Exp7。
- Exp7 在正文中压缩为两行关键结果，完整的 $C_A$、$C_P$、残量和收敛因子保留在原始
  输出中。
- 方法名称统一为“自适应有限步规则”（adaptive finite-step rule）；窗口比较定义为
  “窗口受限离线参考（sampled oracle）”。
- 理论仅以 $m=O(h^{-1})$ 作为固定能量衰减的充分步数尺度；具体比例和尺度比阈值由
  设计组校准。
- 研究报告和理论分析中的说明性语句采用直接陈述，否定表述仅保留在必要的数学结论中。

## 2. 源码复审

- C++ 源码和实验程序共 4315 行；C/C++ 普通注释为 0，shell 脚本除 shebang 外普通
  注释为 0。
- 第二轮逐接口复审确认现存分支均对应算法路径、并行规模选择或可发生的数值故障。
- 8 个运行时异常出口覆盖稀疏 Cholesky 符号分解或 SPD 失败、PCG 正曲率或收敛失败，
  以及结果文件写入失败。
- 七个入口通过 g++ `-fanalyzer`、`-Wall -Wextra -Wpedantic -Werror` 静态分析。

## 3. 构建与完整复跑

- 系统：Linux 6.18.35 x86_64。
- 编译器：g++ 13.3.0，C++17，`-O3 -DNDEBUG`。
- 严格诊断额外包含 shadow、conversion、duplicated-condition、logical-op、
  null-dereference 和 format 检查。
- 当前环境使用 OpenMP 与 pthread，正式运行设为 4 线程。
- 验证脚本通过严格 C++17 直接构建路径编译七个程序。

完整复跑命令：

```bash
TGI_BUILD_DIR=build/final-full TGI_THREADS=4 \
TGI_STEP_TIMEOUT_SECONDS=1800 ./scripts/run_validation.sh full
```

七个程序均正常结束，正式结果均写入 `code/results`，版本字段均为 7.6.0。

| 主题 | 完整复跑核验 |
|---|---|
| 两网格性能与有限路径 | adaptive/global-reference 13/13 收敛，循环和 4459/28972；三条路径复现有限步最小值 |
| 规则质量、稳健性与代价 | gap、seed/RHS 和停止消融均复现；total 中位数 753.910/12152.018 ms |
| 固定物理系数场加密 | 三层共享节点失配 0；adaptive 循环数 122/119/145 |
| 多层初步验证 | adaptive 101/237 次，global-reference 321/770 次 |

## 4. 文档与静态检查

- `theory.tex` 含 80 个唯一标签；全部 `ref/eqref` 和 `cite` 均有对应目标，13 个
  `bibitem` 均被正文引用。
- 在临时 `ctexart` 结构验证类下使用 XeLaTeX 连续编译两遍成功，无未定义引用、悬空
  文献或重复标签。正式排版应在完整中文 TeX 发行版中使用原生 `ctexart`。
- `bash -n code/scripts/run_validation.sh` 通过。
- 全项目检索未发现临时工作标记、旧版本号、旧定位术语或行尾空白。

## 5. 交付核验

交付前执行独立解压、严格直接构建和 quick 运行；压缩包完整性由 `unzip -t` 与 SHA-256
校验确认。
