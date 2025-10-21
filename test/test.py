import sys

sys.path.append(r"/home/zwxdogs/work/USP/py/probe")
sys.path.append(r"/home/zwxdogs/work/USP/py/scatter")

from straight_flow_2d import Straight_flow_2d

scatter_1 = Straight_flow_2d(
    amplitude=1.0,
    N_sca=5,
    N_frame=5,
    velocity=0.01,
    t_step=1e-4,
    width=0.01,
    length=0.1,
)

pos = scatter_1.get_position()

print(pos[3, :, 0])
