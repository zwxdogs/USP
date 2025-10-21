# 超声散射体基类scatter


class Scatter(object):
    def __init__(self, position, amplitude, N_sca, N_frame):
        self.__position = position  # 散射体位置,维度(N_sca, 3, N_frame)
        self.__amplitude = amplitude  # 散射体散射系数
        self.__N_sca = N_sca  # 散射体总数
        self.__N_frame = N_frame  # 总帧数

    def get_position(self):
        return self.__position

    def get_amplitude(self):
        return self.__amplitude

    def get_N_sca(self):
        return self.__N_sca

    def get_N_frame(self):
        return self.__N_frame
