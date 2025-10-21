# 线性阵列
# 集成自探头基类probe
# x坐标沿阵列方向，y坐标垂直于阵列方向，z坐标为阵列面(x-y面)的法向方向
# 坐标原点定义在阵列中心
from probe import Probe


class Linear_probe(Probe):
    def __init__(self, fc, bandwidth, focus, kerf, pitch, N_el, width, height):
        super(Linear_probe, self).__init__(fc, bandwidth, focus)
        self.__kerf = kerf
        self.__pitch = pitch
        self.__N_el = N_el
        self.__width = width  # x方向
        self.__height = height  # y方向

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
