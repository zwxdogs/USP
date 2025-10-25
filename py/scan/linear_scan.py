# Linear scan region.
import scan
from scan.scan_base import Scan_base
import numpy as np
import scan.tools.linear_mesh as linear_mesh


class Linear_scan(Scan_base):
    def __init__(self, min_x, max_x, N_x, min_y, max_y, N_y, min_z, max_z, N_z):
        scan_xyz = linear_mesh.create_mesh(
            min_x, max_x, N_x, min_y, max_y, N_y, min_z, max_z, N_z
        )

        super(Linear_scan, self).__init__(scan_xyz)
