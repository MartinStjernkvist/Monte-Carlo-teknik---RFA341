import scipy.io
from matplotlib.widgets import Slider
from matplotlib.gridspec import GridSpec
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import time
import pandas as pd
from scipy.interpolate import interp1d
from numpy import random

#   ----------------------------------------------------------------------
#   KONSTANTER
#   ----------------------------------------------------------------------

pi = np.pi
E_e = 0.511 * 10 ** 6  # eV
r_e = np.sqrt(0.07941)  # sqrt(b): re2 = e4/Ee2 ≈ 0.07941 b, https://en.wikipedia.org/wiki/Gamma_ray_cross_section
a_0 = 5.29177210903 * 10 ** (-11) * 10 ** (14)  # sqrt(b), bohr radius of hydrogen
c = 3 * 10 ** 8

#   ----------------------------------------------------------------------
#   FUNKTIONER
#   ----------------------------------------------------------------------

def plot_stuff(x_data, y_data, scatter, label_data,
               marker='o', color='blue', x_label='x-label', y_label='y-label', title='1',
               fig_size=(10, 10), symbol_size=100, font_size=30, alpha=1, line_width=10, x_lim=(0, 0), y_lim=(0, 0),
               grid=False, x_scale='linear', y_scale='linear'):
    fig = plt.figure(figsize=fig_size)

    if x_lim != (0, 0) and y_lim != (0, 0):
        plt.ylim(y_lim)
        plt.xlim(x_lim)

    if grid == True:
        plt.grid()

    if scatter == 0:
        plt.plot(x_data, y_data, color=color, alpha=alpha, linewidth=line_width, label=label_data)
    elif scatter == 1:
        plt.scatter(x_data, y_data, marker=marker, color=color, alpha=alpha, s=symbol_size, label=label_data)
    elif scatter == 2:
        plt.plot(x_data, y_data, color=color, alpha=alpha, linewidth=line_width, label=label_data)
        plt.scatter(x_data, y_data, marker=marker, color=color, alpha=alpha, s=symbol_size, label=label_data)
    else:
        for i in range(len(x_data)):
            if scatter[i] == 1:
                plt.scatter(x_data[i], y_data[i], marker=marker[i], color=color[i],
                            alpha=alpha, s=symbol_size, label=label_data[i])
            else:
                plt.plot(x_data[i], y_data[i], color=color[i], alpha=alpha, linewidth=line_width, label=label_data[i])

    font_size_ticks = font_size * 0.85

    plt.xscale(x_scale)
    plt.yscale(y_scale)

    plt.xlabel(x_label, fontsize=font_size)
    plt.ylabel(y_label, fontsize=font_size)

    plt.xticks(fontsize=font_size_ticks)
    plt.yticks(fontsize=font_size_ticks)

    plt.title(title, fontsize=font_size)
    plt.legend(fontsize=font_size_ticks)

    return fig


"""
def scatter_3d(array,
               color='blue', x_label='x-label', y_label='y-label', z_label='z-label', title='1',
               fig_size=(10, 10), symbol_size=100, font_size=30, alpha=1, x_lim=(0, 0), y_lim=(0, 0), z_lim=(0, 0),
               grid=False, x_scale='linear', y_scale='linear', z_scale='linear'):

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(projection='3d')

    # x = np.arange(array.shape[1])
    # y = np.arange(array.shape[0])
    # X, Y = np.meshgrid(x, y)

    x, y, z = array>[:, 0], array[:, 1], array[:, 2]
    ax.scatter(x,y,z, color=color, s = symbol_size, alpha = alpha)

    if x_lim != (0, 0) and y_lim != (0, 0) and z_lim != (0, 0):
        plt.xlim(x_lim)
        plt.ylim(y_lim)
        plt.zlim(z_lim)

    if grid == True:
        plt.grid()

    font_size_ticks = font_size * 0.85

    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_zscale(z_scale)

    ax.set_xlabel(x_label, fontsize=font_size)
    ax.set_ylabel(y_label, fontsize=font_size)
    ax.set_zlabel(z_label, fontsize=font_size)

    ax.set_xticks(fontsize=font_size_ticks)
    ax.set_yticks(fontsize=font_size_ticks)
    ax.set_zticks(fontsize=font_size_ticks)

    ax.set_title(title, fontsize=font_size)
    ax.set_legend(fontsize=font_size_ticks)
    plt.show()
"""