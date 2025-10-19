// 所有散射体的基类scatter
#include <eigen3/Eigen/Dense>

class scatter
{
private:
    int N_sca;
    int N_frame;
    // double amp[];

public:
    scatter(int N_sca, int N_frame)
    {
        this->N_frame = N_frame;
        this->N_sca = N_sca;
    };

    ~scatter() { };
};