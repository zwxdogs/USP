# Class for storage of simulation parameters


class Simu_para(object):
    def __init__(self, c0, fs):
        self.__c0 = c0  # speed of sound
        self.__fs = fs  # sampling frequency
        self.__dt = 1 / fs  # time step

    def get_c0(self):
        return self.__c0

    def get_fs(self):
        return self.__fs

    def get_dt(self):
        return self.__dt
