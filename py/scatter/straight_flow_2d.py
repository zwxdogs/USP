# 二维沿直线流动散射体分布
# 坐标原点定义在入口中心
from scipy import rand
from scatter import Scatter
import numpy as np
import random


class Straight_flow_2d(Scatter):
    def calc_position(self, N_sca, N_frame, velocity, t_step, width, length):
        # 初始坐标计算
        position = np.zeros((N_sca, 3, N_frame))
        y0 = np.zeros(N_sca)
        position[:, 1, 0] = y0
        for n in range(N_sca):
            random_x = random.random()
            random_z = random.uniform(-0.5, 0.5)
            position[n, 0, 0] = random_x * length
            position[n, 2, 0] = random_z * width

        # 各时间帧坐标计算
        for f in range(1, N_frame):
            for n in range(N_sca):
                position[n, 0, f] = position[n, 0, f - 1] + velocity * t_step
                # 周期进出口
                if position[n, 0, f] > length:
                    position[n, 0, f] = position[n, 0, f] - length

                position[n, 1, f] = position[n, 1, 0]  # y坐标不变
                position[n, 2, f] = position[n, 2, 0]  # z坐标不变

        return position

    def __init__(self, amplitude, N_sca, N_frame, velocity, t_step, width, length):
        self.__velocity = velocity  # 散射体流动速度
        self.__t_step = t_step  # 时间步长
        self.__width = width  # 流动区域宽度
        self.__length = length  # 流动区域长度

        position = self.calc_position(N_sca, N_frame, velocity, t_step, width, length)

        super().__init__(position, amplitude, N_sca, N_frame)

    def get_velocity(self):
        return self.__velocity

    def get_t_step(self):

        return self.__t_step

    def get_width(self):
        return self.__width

    def get_length(self):
        return self.__length
