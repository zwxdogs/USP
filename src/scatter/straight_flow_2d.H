// 沿直线流动的散射体分布

#include "scatter.H"

class straight_flow_2d : public scatter
{
private:
    vec_d amp;
    mat_d pos;
    double v;
    double t_step;
    double width;

public:
    // Constructor
    straight_flow_2d(int N_sca, int N_frame, double v, double t_step, double width) :
        scatter(N_sca, N_frame)
    {
        this->v = v;
        this->t_step = t_step;
        this->width = width;
    }

    // Destructor
    ~straight_flow_2d() { };
};