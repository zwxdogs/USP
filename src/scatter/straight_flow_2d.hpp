// 沿直线流动的散射体分布

#include "scatter.hpp"
#include <random>

class straight_flow_2d : public scatter
{
private:
    double v;
    double t_step;
    double width;
    double length;

    // Member functions

    // 计算散射体的位置
    void calc_pos()
    {
        // 随机数生成器
        std::random_device seed;
        std::ranlux48 engine(seed());
        std::uniform_real_distribution<double> distrib_1(-0.5, 0.5);
        std::uniform_real_distribution<double> distrib_2(0, 1);
        // 初始坐标计算
        vec_d y0 = vec_d::Zero(this->N_sca);
        vec_d x0(this->N_sca);
        vec_d z0(this->N_sca);
        for ( int i = 0; i < this->N_sca; i++ ) {
            double random_x = distrib_2(engine);
            double random_z = distrib_1(engine);
            x0(i) = random_x * length;
            y0(i) = random_z * width;
        }
        // 各时间帧位置计算
        for ( int i = 0; i < this->N_frame; i++ ) {
            vec_d tmp_amp = vec_d::Ones(this->N_sca);  // 暂时只需要散射系数为1。
            mat_d tmp_pos = mat_d::Zero(this->N_sca, 3);
            for ( int j = 0; j < this->N_sca; j++ ) {
                double y = y0(j);
                double x = x0(j) + this->v * this->t_step * i;  // 以固定速度沿x方向流动。
                if ( x > length ) {
                    x = x - length;  // 周期边界。
                }
                double z = z0(j);  // z方向不变。
                tmp_pos(j, 0) = x;
                tmp_pos(j, 1) = y;
                tmp_pos(j, 2) = z;
            }
            this->amp.emplace_back(tmp_amp);
            this->pos.emplace_back(tmp_pos);
        }
    }

public:
    // Constructor
    straight_flow_2d(int N_sca, int N_frame, double v, double t_step, double width, double length) :
        scatter(N_sca, N_frame)
    {
        this->v = v;
        this->t_step = t_step;
        this->width = width;
        this->length = length;

        calc_pos();
    }

    // Destructor
    ~straight_flow_2d() { };

    // Member functions
};