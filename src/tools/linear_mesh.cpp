// Calculate linear mesh (x, y, z) by three one-dimensional arrays.
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <eigen3/Eigen/Dense>

namespace py = pybind11;
using matXd = Eigen::MatrixXd;

py::array_t<double>
value_double(py::array_t<double>& arr1)
{
    auto r1 = arr1.unchecked<2>();

    py::array_t<double> out = py::array_t<double>(arr1.size());
    out.resize({arr1.shape()[0], arr1.shape()[1]});
    auto r_out = out.mutable_unchecked<2>();

    for ( int i = 0; i < arr1.shape()[0]; i++ ) {
        for ( int j = 0; j < arr1.shape()[1]; j++ ) {
            r_out(i, j) = r1(i, j) * 2.0;
        }
    }

    return out;
}

PYBIND11_MODULE(test_np, m)
{
    m.doc() = "Numpy array value double function module";

    m.def("value_double", &value_double, "A function that doubles the input array values");
}