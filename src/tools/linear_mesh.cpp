// Calculate linear mesh grid.

#include <nanobind/nanobind.h>  // nanobind模块
#include <nanobind/ndarray.h>   // 数组模块

namespace nb = nanobind;

// using arr_3c = nb::ndarray<double, nb::numpy, nb::c_contig, nb::shape<-1, 3>>;

nb::ndarray<double, nb::numpy, nb::c_contig>
create_mesh(double min_x, double max_x, int N_x, double min_y, double max_y, int N_y, double min_z, double max_z, int N_z)
{
    double dx = (max_x - min_x) / (N_x - 1);
    double dy = (max_y - min_y) / (N_y - 1);
    double dz = (max_z - min_z) / (N_z - 1);

    double* data = new double[N_x * N_y * N_z * 3];

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

    return nb::ndarray<double, nb::numpy, nb::c_contig>(
      data, shape, owner);
}

NB_MODULE(linear_mesh, m)
{
    m.doc() = "Calculate linear mesh grid.";

    m.def("create_mesh", &create_mesh, "Calculate linear mesh grid.", nb::rv_policy::take_ownership);  // 绑定函数
}