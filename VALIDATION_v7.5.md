# multiscale_two_grids_iteration v7.5 验证记录

数值源码版本：7.5.0

理论与研究报告版本：v7.5

## 1. 精简与一致性

- C++ 源码和实验程序由 4707 行精简为 4315 行，共减少 392 行（8.3%）。
- C/C++ 普通注释为 0；shell 脚本除 shebang 外普通注释为 0。
- 删除未使用的插值报告结构、列迭代统计、冗余矩阵接口、重复维数与层次检查，以及由
  内部构造保证不会触发的分支。
- 保留 8 个运行时异常出口，分别覆盖稀疏 Cholesky 符号分解或 SPD 失败、PCG 正曲率或
  收敛失败，以及结果文件写入失败。
- 项目术语统一为“自适应有限步规则”（adaptive finite-step rule）。理论给出
  $m=O(h^{-1})$ 的充分步数尺度，规则中的比例和阈值由设计组校准。

## 2. 构建与运行环境

- 系统：Linux 6.18.35 x86_64。
- 编译器：g++ 13.3.0，C++17，`-O3 -DNDEBUG`。
- 诊断：`-Wall -Wextra -Wpedantic -Werror`，并启用编译器支持的 shadow、conversion、
  duplicated-condition、logical-op、null-dereference 和 format 诊断。
- 并行：OpenMP 与 pthread，正式运行使用 4 线程。
- 当前环境未提供 CMake；验证脚本使用等价的严格 C++17 直接构建路径完成七个程序。

## 3. 完整数值复跑

命令：

```bash
TGI_BUILD_DIR=build/final-full TGI_THREADS=4 \
TGI_STEP_TIMEOUT_SECONDS=1800 ./scripts/run_validation.sh full
```

七个程序均正常结束，七份正式结果均写入 `code/results`，版本字段均为 7.5.0。

| 实验 | 完整复跑核验 |
|---|---|
| 两网格主比较 | adaptive/global-reference 13/13 收敛，累计循环数 4459/28972 |
| 有限 PCG 路径 | 三条路径均复现能量下降与两网格循环数的有限步最小值 |
| sampled oracle | 设计组平均/最大 gap 8.93%/29.65%，验证组 27.34%/54.20% |
| seed、RHS 与计时 | seed 循环和 327/866，RHS 循环和 469/1507；total 中位数 570.872/8901.375 ms |
| 停止策略消融 | 收敛数 6/6、5/6、5/6；记录循环和 1147/21098/25658 |
| 固定物理加密 | 三层共享节点失配数 0；adaptive 循环数 122/119/145 |
| 三层 V-cycle | adaptive 101/237 次，global-reference 321/770 次 |

## 4. 文档与静态检查

- `theory.tex` 含 80 个唯一标签；全部 `ref/eqref` 和 `cite` 均有对应目标，13 个
  `bibitem` 均被正文引用。
- 在临时 `ctexart` 结构验证类下使用 XeLaTeX 连续编译两遍成功，无未定义引用、悬空
  文献或重复标签。正式排版仍应在完整中文 TeX 发行版中使用原生 `ctexart`。
- `bash -n code/scripts/run_validation.sh` 通过。
- 全项目检索未发现临时工作标记、旧版本号、旧定位术语或行尾空白。

## 5. 交付核验

交付前执行独立解压、严格直接构建和 quick 运行；压缩包完整性由 `unzip -t` 与 SHA-256
校验确认。
