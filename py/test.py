# 测试模拟
import time
from probe import Rect_probe
from scan import Linear_scan
from simulation import Simu_para
from simulation import calc_SIR
from plot import plot_SIR

para = Simu_para(1600, 100e6)
prob = Rect_probe(5e6, 0.8, [0, 0, 0], 3e-3, 5e-3)
prob.set_subapertures(10, 10)

sca = Linear_scan(0e-3, 6e-3, 4, 0e-3, 0e-3, 1, 0e-3, 5e-3, 2)

simulation = calc_SIR(prob, sca, para)

plot_SIR(simulation, para, 5)
