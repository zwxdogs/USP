// Calculate linear mesh grid.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;
using arr_np_c = nb::ndarray<double, nb::numpy, nb::c_contig>;

arr_np_c
create_mesh(double min_x, double max_x, int N_x, double min_y, double max_y, int N_y, double min_z, double max_z, int N_z)
{
    if ( N_x <= 0 || N_y <= 0 || N_z <= 0 ) {
        throw std::runtime_error("Number of points in each direction must be positive.");
    }
    double dx = 0;
    double dy = 0;
    double dz = 0;
    if ( N_x != 1 ) {
        dx = (max_x - min_x) / (N_x - 1);  // One point handling.
    }
    if ( N_y != 1 ) {
        dy = (max_y - min_y) / (N_y - 1);
    }
    if ( N_z != 1 ) {
        dz = (max_z - min_z) / (N_z - 1);
    }

    double* data = new double[N_x * N_y * N_z * 3];
    // create mesh grid
    size_t index = 0;
    for ( int i = 0; i < N_x; ++i ) {
        double x_now = min_x + i * dx;
        for ( int j = 0; j < N_y; ++j ) {
            double y_now = min_y + j * dy;
            for ( int k = 0; k < N_z; ++k ) {
                data[index++] = x_now;
                data[index++] = y_now;
                data[index++] = min_z + k * dz;
            }
        }
    }

    nb::capsule owner(data,
                      [](void* p) noexcept {  // 析构函数
                          delete[] (double*)p;
                      });

    std::initializer_list<size_t> shape = {(size_t)(N_x * N_y * N_z), 3};

    return arr_np_c(
      data, shape, owner);
}

NB_MODULE(linear_mesh, m)
{
    m.doc() = "Calculate linear mesh grid.";

    m.def("create_mesh", &create_mesh, "Calculate linear mesh grid.", nb::rv_policy::take_ownership);  // 绑定函数
}