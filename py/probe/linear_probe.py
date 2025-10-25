# Linear array probe class.
# x coordinates along the array direction, y coordinates perpendicular to the array direction, z coordinates normal to the array plane (x-y plane).
# Origin is at the center of the probe.
from turtle import position
from probe.probe_base import Probe
import numpy as np


class Linear_probe(Probe):
    def calc_position(pitch, N_el):
        position = np.zeros((N_el, 3))
        for n in range(N_el):
            position[n, 0] = n * pitch - (N_el - 1) * pitch / 2  # x coordinate
            position[n, 1] = 0  # y coordinate
            position[n, 2] = 0  # z coordinate

        return position

    def __init__(self, fc, bandwidth, focus, kerf, pitch, N_el, height):
        position = Linear_probe.calc_position(pitch, N_el)
        super(Linear_probe, self).__init__(fc, bandwidth, focus, position)
        self.__kerf = kerf
        self.__pitch = pitch
        self.__N_el = N_el
        self.__width = pitch - kerf
        self.__height = height
        self.__symbol = self.symbol_list[1]  # "linear"

    def get_kerf(self):
        return self.__kerf

    def get_pitch(self):
        return self.__pitch

    def get_N_el(self):

        return self.__N_el

    def get_width(self):
        return self.__width

    def get_height(self):
        return self.__height

    def get_symbol(self):
        return self.__symbol
