# 所有探头的基类probe

import numpy as np


class Probe(object):
    def __init__(self, fc, bandwidth, focus, center=np.array([0, 0, 0])):
        self.__fc = fc
        self.__bandwidth = bandwidth
        self.__focus = focus  # 数组，包含x y z三个坐标。
        self.__center = center  # 数组，包含x y z三个坐标，默认为(0, 0, 0)

    def get_fc(self):
        return self.__fc

    def get_bandwidth(self):
        return self.__bandwidth

    def get_focus(self):
        return self.__focus

    def get_center(self):
        return self.__center
