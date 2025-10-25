// Calculate Spatial Impulse Response (SIR).

#include <nanobind/nanobind.h>  // nanobind模块
#include <nanobind/ndarray.h>   // 数组模块

namespace nb = nanobind;
using arr_np_c = nb::ndarray<double, nb::numpy, nb::c_contig>;

arr_np_c
calc(arr_np_c probe_line,
     arr_np_c probe_corners,
     arr_np_c scan_position,
     double dt,
     double c0)
{
    // To be implemented
    // 先找最近点和最远点，设定t的范围。
    // 在t范围内循环，确定交点->确定弧段->加入总和
    return arr_np_c();
}

NB_MODULE(SIR_calc, m)
{
    m.doc() = "Calculate Spatial Impulse Response (SIR).";

    m.def("calc", &calc, "Calculate Spatial Impulse Response (SIR).", nb::rv_policy::take_ownership);  // 绑定函数
}