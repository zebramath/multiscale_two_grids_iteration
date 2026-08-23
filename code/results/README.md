# Results

`experiment1`--`experiment6` 是随 v2.11.0 保留的统一基线结果，因此文本中的
版本号仍为 2.11.0；其源码已随 v3.0.0 一起编译通过，重新运行后会写出 3.0.0
版本标识。

`experiment7`--`experiment9` 是 v3.0.0 的新增实测结果：

- experiment7：128×128 原始通道上的自适应 PCG；
- experiment8：四类系数场上的新旧支撑策略消融；
- experiment9：18 问题完整鲁棒性矩阵。

CSV 用于后续统计，TXT 为相同结果的可读报告。墙钟时间来自单次运行，不应跨机器
比较绝对值。
