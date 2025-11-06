# 测试模拟
import time
from probe import Rect_probe
from scan import Linear_scan
from simulation import Simu_para
from simulation import calc_SIR
import plot
from simulation.tools.apodization import apo_types

para = Simu_para(1600, 100e6)
prob = Rect_probe(5e6, 1, [0, 0, 0], 5e-3, 5e-3)
prob.set_subapertures(10, 10)

sca = Linear_scan(0e-3, 5e-3, 5, 0e-3, 0e-3, 1, 5e-3, 5e-3, 1)

simulation = calc_SIR(prob, sca, para, apo_types.NONE)

plot.plot_SIR(simulation, para, 0)
# plot.plot_SIR(simulation, para, 1)
# plot.plot_SIR(simulation, para, 2)
# plot.plot_SIR(simulation, para, 3)
# plot.plot_SIR(simulation, para, 4)
