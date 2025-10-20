// 所有散射体的基类scatter
#include <eigen3/Eigen/Dense>
#include "../tools/VTK/point_cloud.hpp"

using vec_d = Eigen::VectorXd;
using vec_v = std::vector<vec_d>;
using mat_d = Eigen::MatrixXd;
using vec_m = std::vector<mat_d>;

class scatter
{
protected:
    vec_v amp;
    vec_m pos;
    int N_sca;
    int N_frame;

public:
    // Constructor
    scatter(int N_sca, int N_frame)
    {
        this->N_frame = N_frame;
        this->N_sca = N_sca;
    };

    // Destructor
    ~scatter() { };

    // Member functions
    vec_m get_pos()
    {
        return this->pos;
    };

    // 绘制点云
    void plot_pos(int frame) {
        // 通过pybind11将点数据输出到python，再通过matplotlib画图。

    };
};