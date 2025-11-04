# Use matplotlib to plot SIR results.
import matplotlib.pyplot as plt
import numpy as np


def plot_SIR(SIR_result, para, sca_index):
    max_step = len(SIR_result[sca_index])
    time_axis = np.linspace(0, (max_step - 1) * para.get_dt(), max_step)

    plt.plot(time_axis, SIR_result[sca_index])
    plt.show()
