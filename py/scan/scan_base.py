# Base class of all of the scan types.


class Scan_base:
    def __init__(self, scan_xyz):
        self.__scan_xyz = scan_xyz

    def get_scan_xyz(self):
        return self.__scan_xyz
