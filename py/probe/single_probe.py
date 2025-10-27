# Single probe class.
from probe.probe_base import Probe
import numpy as np


class Single_probe(Probe):
    def __init__(self, fc, bandwidth, focus, width, height):
        super(Single_probe, self).__init__(fc, bandwidth, focus, [0, 0, 0])
        self.__symbol = self.symbol_list[0]  # "single"
        self.__width = width  # x direction
        self.__height = height  # y direction
        self.__N_el = 1
        # Polygon aperture
        self.__line = np.zeros((4, 3))
        self.__line = np.array(
            [
                [0, 1, -height / 2],  # line_1: y = height/2
                [1, 0, -width / 2],  # line_2: x = width/2
                [0, 1, height / 2],  # line_3: y = -height/2
                [1, 0, width / 2],  # line_4: x = -width/2
            ]
        )
        self.__corners = [
            [width / 2, height / 2, 0],  # line1 & line2 intersection
            [width / 2, -height / 2, 0],  # line2 & line3 intersection
            [-width / 2, -height / 2, 0],  # line3 & line4 intersection
            [-width / 2, height / 2, 0],  # line4 & line1 intersection
        ]

    def get_symbol(self):
        return self.__symbol

    def get_width(self):
        return self.__width

    def get_height(self):
        return self.__height

    def get_line(self):
        return self.__line

    def get_corners(self):
        return self.__corners

    def get_N_el(self):
        return self.__N_el
