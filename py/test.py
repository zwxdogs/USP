# 测试模拟
import time
from probe import Rect_probe
from scan import Linear_scan
from simulation import Simu_para
from simulation import calc_SIR
import matplotlib.pyplot as plt
import numpy as np

para = Simu_para(1540, 100e6)
prob = Rect_probe(5e6, 0.8, [0, 0, 0], 4e-3, 5e-3)
prob.set_subapertures(10, 10)

sca = Linear_scan(-3e-3, 3e-3, 7, -2e-3, 2e-3, 5, 5e-3, 6e-3, 2)

simulation = calc_SIR(prob, sca, para)
max_step = len(simulation[0])
time_axis = np.linspace(0, (max_step - 1) * para.get_dt(), max_step)

plt.plot(time_axis, simulation[0])
plt.show()
