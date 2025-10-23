# Base class for all of the ultrasound scatterers.
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Scatter(object):
    def __init__(self, position, amplitude, N_sca, N_frame):
        self._position = position  # dimension (N_sca, 3, N_frame)
        self._amplitude = amplitude  # scatterer amplitude
        self._N_sca = N_sca
        self._N_frame = N_frame

    # View scatterer distribution in 3D by matplotlib
    def plot_scatter(self, frame_idx):
        fig = plt.figure()
        # 3d coordinate system
        ax = fig.add_subplot(111, projection="3d")
        # 3d scatter plot
        ax.scatter(
            self._position[:, 0, frame_idx],
            self._position[:, 1, frame_idx],
            self._position[:, 2, frame_idx],
            c="r",
            marker="o",
        )
        ax.set_xlabel("X Label")
        ax.set_ylabel("Y Label")
        ax.set_zlabel("Z Label")
        plt.show()

    def get_position(self):
        return self._position

    def get_amplitude(self):
        return self._amplitude

    def get_N_sca(self):
        return self._N_sca

    def get_N_frame(self):
        return self._N_frame
