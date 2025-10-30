import SIR_calc

# 测试模拟
from probe import Rect_probe
from scan import Linear_scan
from simulation import Simu_para
import matplotlib.pyplot as plt
import numpy as np

para = Simu_para(1540, 100e6)
prob = Rect_probe(5e6, 0.8, [0, 0, 0], 2e-3, 3e-3)
prob.set_subapertures(2, 3)

sca = Linear_scan(0e-3, 0e-3, 1, 0, 0, 1, 5e-3, 5e-3, 1)

simulation = SIR_calc.calc_polygon(
    prob.get_line(),
    prob.get_calc_line(),
    prob.get_corners(),
    prob.get_calc_corners(),
    sca.get_scan_xyz(),
    para.get_dt(),
    para.get_c0(),
)

# print("Simulation Result:")
# print(simulation[0])
# print("\nSize of result:")
# print(len(simulation[0]))


plt.plot(simulation[0])
plt.show()
