# Base class of all of the probes.

import numpy as np
from abc import abstractmethod


class Probe(object):
    def __init__(self, fc, bandwidth, focus, position, center=np.array([0, 0, 0])):
        self._fc = fc
        self._bandwidth = bandwidth
        self._focus = focus  # array (x, y, z)
        self._center = center  # array (x, y, z), default at origin
        self.symbol_list = ["rect", "linear"]
        self._position = position  # array (x, y, z)

    def get_fc(self):
        return self._fc

    def get_bandwidth(self):
        return self._bandwidth

    def get_focus(self):
        return self._focus

    def get_center(self):
        return self._center

    def get_symbol_list(self):
        return self.symbol_list

    def get_position(self):
        return self._position

    @abstractmethod
    def get_symbol(self):
        pass
