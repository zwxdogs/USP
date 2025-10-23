# 2d straight flow scatterer class.
# Origin is at the center of the inlet.
# Noted that it's need to transform the coordinates to adapt probes
# when used in simulation.
from scatter.scatter_base import Scatter
import numpy as np
import random


class Straight_flow_2d(Scatter):
    def calc_position(self, N_sca, N_frame, velocity, t_step, width, length):
        # time frame 0 position
        position = np.zeros((N_sca, 3, N_frame))
        y0 = np.zeros(N_sca)
        position[:, 1, 0] = y0
        for n in range(N_sca):
            random_x = random.random()
            random_z = random.uniform(-0.5, 0.5)
            position[n, 0, 0] = random_x * length
            position[n, 2, 0] = random_z * width

        # Calculate frames position at all the time steps
        for f in range(1, N_frame):
            for n in range(N_sca):
                position[n, 0, f] = position[n, 0, f - 1] + velocity * t_step
                # Cycle boundary
                if position[n, 0, f] > length:
                    position[n, 0, f] = position[n, 0, f] - length

                position[n, 1, f] = position[n, 1, 0]  # y coordinate unchanged
                position[n, 2, f] = position[n, 2, 0]  # z coordinate unchanged

        return position

    def __init__(self, amplitude, N_sca, N_frame, velocity, t_step, width, length):
        self.__velocity = velocity
        self.__t_step = t_step  # it's related to the effctive PRF.
        self.__width = width  # flow region width
        self.__length = length  # flow region length

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
