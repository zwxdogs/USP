# USP
目标是构建一个超声代码测试流程，从而为实验算法进行前期测试。
- 采用Python开发。
- 采用C++进行复杂计算。
- 通过nanobind进行绑定。

## 记录
- 散射体坐标采用Numpy矩阵保存。
- 使用空间脉冲响应方法。
- SIR计算结果使用`std::vector<std::vector<double>>`保存。
 - 返回到Python为List嵌套。
