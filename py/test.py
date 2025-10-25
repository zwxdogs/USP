# from scan import Linear_scan

# sca = Linear_scan(-1, 1, 5, -2, 2, 5, -3, 3, 5)
# print(sca.get_scan_xyz())

import numpy as n
from probe import Linear_probe

lp = Linear_probe(5e6, 0.6, [0, 0, 0], 0.3e-3, 0.4e-3, 5, 5e-3)

print(lp.get_position())
