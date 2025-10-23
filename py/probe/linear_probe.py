# Linear array probe class.
# x coordinates along the array direction, y coordinates perpendicular to the array direction, z coordinates normal to the array plane (x-y plane).
# Origin is at the center of the probe.
from probe.probe_base import Probe


class Linear_probe(Probe):
    def __init__(self, fc, bandwidth, focus, kerf, pitch, N_el, width, height):
        super(Linear_probe, self).__init__(fc, bandwidth, focus)
        self.__kerf = kerf
        self.__pitch = pitch
        self.__N_el = N_el
        self.__width = width
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
