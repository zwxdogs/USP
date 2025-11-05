# Single Rectangle probe class.

# line store: y - kx + b and x = c. use array[N_lines, 3]. First column: 0 or 1 for define whether the slope is infinte;
# second column: k or 0, depending on slope; third column: b or c, depending on slope.

from probe.probe_base import Probe
import numpy as np


class Rect_probe(Probe):
    def __init__(self, fc, bandwidth, focus, width, height):
        super(Rect_probe, self).__init__(fc, bandwidth, focus, [0, 0, 0])
        self.__symbol = self.symbol_list[0]  # "single"
        self.__width = width  # x direction
        self.__height = height  # y direction
        self.__N_el = 1
        # Polygon aperture
        self.__line = np.zeros((4, 3))
        self.__line = np.array(
            [
                [0, 0, height / 2],  # line_1: y = height/2
                [1, 0, width / 2],  # line_2: x = width/2
                [0, 0, -height / 2],  # line_3: y = -height/2
                [1, 0, -width / 2],  # line_4: x = -width/2
            ]
        )
        self.__corners = np.array(
            [
                [width / 2, height / 2, 0],  # line1 & line2 intersection
                [width / 2, -height / 2, 0],  # line2 & line3 intersection
                [-width / 2, -height / 2, 0],  # line3 & line4 intersection
                [-width / 2, height / 2, 0],  # line4 & line1 intersection
            ]
        )

    # Divide subapertures
    def set_subapertures(self, N_sub_x, N_sub_y):
        self.__N_sub_x = N_sub_x
        self.__N_sub_y = N_sub_y
        self.__inside_line = np.zeros(((N_sub_x - 1) + (N_sub_y - 1), 3))
        self.__inside_corners = np.zeros(
            ((N_sub_x - 1) * 2 + (N_sub_x + 1) * (N_sub_y - 1), 3)
        )
        for i in range(1, N_sub_x):
            self.__inside_line[i - 1, 0] = 1
            self.__inside_line[i - 1, 1] = 0
            x_position = -self.__width / 2 + self.__width / N_sub_x * i
            self.__inside_line[i - 1, 2] = x_position

            self.__inside_corners[i - 1, 0] = x_position
            self.__inside_corners[i - 1, 1] = self.__height / 2
            self.__inside_corners[i - 1, 2] = 0
            self.__inside_corners[i - 1 + N_sub_x - 1, 0] = x_position
            self.__inside_corners[i - 1 + N_sub_x - 1, 1] = -self.__height / 2
            self.__inside_corners[i - 1 + N_sub_x - 1, 2] = 0

        for i in range(1, N_sub_y):
            self.__inside_line[N_sub_x - 1 + i - 1, 0] = 0
            self.__inside_line[N_sub_x - 1 + i - 1, 1] = 0
            y_position = -self.__height / 2 + self.__height / N_sub_y * i
            self.__inside_line[N_sub_x - 1 + i - 1, 2] = y_position

            for j in range(N_sub_x + 1):
                x_position = -self.__width / 2 + self.__width / N_sub_x * j
                self.__inside_corners[
                    (N_sub_x - 1) * 2 + (N_sub_x + 1) * (i - 1) + j, 0
                ] = x_position
                self.__inside_corners[
                    (N_sub_x - 1) * 2 + (N_sub_x + 1) * (i - 1) + j, 1
                ] = y_position
                self.__inside_corners[
                    (N_sub_x - 1) * 2 + (N_sub_x + 1) * (i - 1) + j, 2
                ] = 0

        self.__calc_line = np.concatenate((self.__line, self.__inside_line), axis=0)
        self.__calc_corners = np.concatenate(
            (self.__corners, self.__inside_corners), axis=0
        )

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

    def get_N_sub_x(self):
        return self.__N_sub_x

    def get_N_sub_y(self):
        return self.__N_sub_y

    def get_inside_line(self):
        return self.__inside_line

    def get_inside_corners(self):
        return self.__inside_corners

    def get_calc_line(self):
        return self.__calc_line

    def get_calc_corners(self):
        return self.__calc_corners
