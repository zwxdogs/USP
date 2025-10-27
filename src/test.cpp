#include <vector>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>  // 关键头文件，用于自动转换std::vector

namespace nb = nanobind;

// 假设这个函数在循环中生成不同长度的向量
std::vector<std::vector<int>>
generate_ragged_data()
{
    std::vector<std::vector<int>> all_data;

    // 模拟一个循环，例如循环5次
    for ( int i = 0; i < 5; ++i ) {
        std::vector<int> current_loop_data;
        int random_length = 2 + (i % 4);  // 产生不同长度：2, 3, 4, 5, 2

        for ( int j = 0; j < random_length; ++j ) {
            current_loop_data.push_back(i * 10 + j);
        }
        all_data.push_back(current_loop_data);
    }
    return all_data;  // 直接返回对象，nanobind会处理转换
}

// 使用nanobind进行绑定
NB_MODULE(test, m)
{
    m.doc() = "A simple example of passing vector of vectors to Python";
    m.def("generate_ragged_data", &generate_ragged_data, "Generates and returns a list of lists with variable lengths.");
}