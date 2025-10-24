# Linear scan region.
from scan.scan_base import Scan_base
import numpy as np


class Linear_scan(Scan_base):
    def __init__(self, min_x, max_x, N_x, min_y, max_y, N_y, min_z, max_z, N_z):
        scan_x = np.linspace(min_x, max_x, N_x).reshape(N_x, 1)
        scan_y = np.linspace(min_y, max_y, N_y).reshape(N_y, 1)
        scan_z = np.linspace(min_z, max_z, N_z).reshape(N_z, 1)

        super(Linear_scan, self).__init__(scan_xyz)
