// Calculate linear mesh (x, y, z) by three one-dimensional arrays.
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <eigen3/Eigen/Dense>

namespace nb = nanobind;

using matXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using test_arr = nb::ndarray<double, nb::numpy, nb::shape<-1, 3>>;

test_arr
value_double(test_arr arr1)
{
    int row = arr1.shape(0);
    int col = arr1.shape(1);

    for ( int i = 0; i < row; i++ ) {
        for ( int j = 0; j < col; j++ ) {
            arr1(i, j) = arr1(i, j) * 2.0;
        }
    }
    return arr1;
}

NB_MODULE(test_np, m)
{
    m.doc() = "Numpy array value double function module";

    m.def("value_double", &value_double, "A function that doubles the input array values");
}