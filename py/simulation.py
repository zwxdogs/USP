# 测试模拟
from probe import Single_probe
from scan import Linear_scan
from simulation import Simu_para

para = Simu_para(1540, 50e6)
prob = Single_probe(5e6, 0.8, [0, 0, 0], 3e-3, 5e-3)
sca = Linear_scan(0, 0, 1, 0, 0, 1, 5e-3, 5e-3, 1)
