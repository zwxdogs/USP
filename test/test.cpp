#include "../src/scatter/straight_flow_2d.hpp"
#include <iostream>

int
main()
{
    int N_sca = 1000;
    int N_frame = 50;
    double v = 0.01;
    double t_step = 0.1;
    double width = 1.0;
    double length = 5.0;

    straight_flow_2d scatterers(N_sca, N_frame, v, t_step, width, length);
    // vec_m positions = scatterers.get_pos();
    // // 输出每一帧的第一个散射体位置作为验证
    // for ( int i = 0; i < N_frame; i++ ) {
    //     std::cout << "Frame " << i << ": ("
    //               << positions[i](0, 0) << ", "
    //               << positions[i](0, 1) << ", "
    //               << positions[i](0, 2) << ")" << std::endl;
    // }

    scatterers.plot_pos(1);

    return 0;
}