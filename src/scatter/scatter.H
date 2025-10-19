// 所有散射体的基类scatter
#include <eigen3/Eigen/Dense>

using vec_d = Eigen::VectorXd;
using mat_d = Eigen::MatrixXd;

class scatter
{
private:
    int N_sca;
    int N_frame;

public:
    scatter(int N_sca, int N_frame)
    {
        this->N_frame = N_frame;
        this->N_sca = N_sca;
    };

    ~scatter() { };
};