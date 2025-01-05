# Genom att skriva "from imports import *" i nya filer kommer alla paket, konstanter och funktioner
# från denna fil att importeras till den nya filen.

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
from scipy.optimize import curve_fit
import multiprocessing as mp
import tkinter as tk
from tkinter import simpledialog
from numba import jit

#   ----------------------------------------------------------------------
#   KONSTANTER
#   ----------------------------------------------------------------------

pi = np.pi
E_e = 0.511 * 10 ** 6  # eV
r_e = np.sqrt(0.07941)  # sqrt(b): re2 = e4/Ee2 ≈ 0.07941 b, https://en.wikipedia.org/wiki/Gamma_ray_cross_section
a_0 = 5.29177210903 * 10 ** (-11) * 10 ** (14)  # sqrt(b), bohr radius of hydrogen
c = 3 * 10 ** 8


# Neutroner???
massa_H = 1 # enhet u
# massa_C = 6
# massa_N = 7
# massa_O = 8
# massa_Na = 11
# massa_Mg = 12
# massa_P = 15
# massa_S = 16
# massa_K = 19
# massa_Ca = 20
massa_C = 12
massa_N = 14
massa_O = 16
massa_Na = 23
massa_Mg = 24.3
massa_P = 31
massa_S = 32
massa_K = 39
massa_Ca = 40

# Sammansättning av vävnad, sida 71 tabeller radiofysik del 1 canvas RFA331.
vävnad_sammansättning_vektor = np.array(
    [10.2, 12.3 / massa_C, 3.5 / massa_N, 72.9 / massa_O, 0.08 / massa_Na, 0.02 / massa_Mg, 0.2 / massa_P,
     0.5 / massa_S, 0.3 / massa_K, 0.007 / massa_Ca])

# Atomic-Electron Binding Energies.pdf canvas RFA331.
K_alpha_H = 0.0136 * 10 ** 3
K_alpha_C = 0.2838 * 10 ** 3
K_alpha_N = 0.4016 * 10 ** 3
K_alpha_O = 0.5320 * 10 ** 3
K_alpha_Na = 1.0721 * 10 ** 3
K_alpha_Mg = 1.3050 * 10 ** 3
K_alpha_P = 2.1455 * 10 ** 3
K_alpha_S = 2.4720 * 10 ** 3
K_alpha_K = 3.6074 * 10 ** 3
K_alpha_Ca = 4.0381 * 10 ** 3

K_alpha_vektor = np.array(
    [K_alpha_H, K_alpha_C, K_alpha_N, K_alpha_O, K_alpha_Na, K_alpha_Mg, K_alpha_P, K_alpha_S, K_alpha_K, K_alpha_Ca])

# Viktat medelvärde
K_alpha = (1 / np.sum(vävnad_sammansättning_vektor)) * vävnad_sammansättning_vektor @ K_alpha_vektor.T
# print(K_alpha)

foton_energi_threshhold = K_alpha_H

# Fluorescence and Coster-Kronig Yields, sida 9 tabeller radiofysik del 1 canvas RFA331.
fluorescence_yield_H = 0
fluorescence_yield_C = 0.0028
fluorescence_yield_N = 0.0052
fluorescence_yield_O = 0.0083
fluorescence_yield_Na = 0.023
fluorescence_yield_Mg = 0.030
fluorescence_yield_P = 0.063
fluorescence_yield_S = 0.078
fluorescence_yield_K = 0.140
fluorescence_yield_Ca = 0.163

fluorescence_yield_vektor = np.array(
    [fluorescence_yield_H, fluorescence_yield_C, fluorescence_yield_N, fluorescence_yield_O, fluorescence_yield_Na,
     fluorescence_yield_Mg, fluorescence_yield_P, fluorescence_yield_S, fluorescence_yield_K, fluorescence_yield_Ca])

# Viktat medelvärde.
fluorescence_yield = (1 / np.sum(
    vävnad_sammansättning_vektor)) * vävnad_sammansättning_vektor @ fluorescence_yield_vektor.T

# Sidlängd på voxlarna i matriserna.
voxel_sidlängd = 0.15  # cm


Lu177_energi = [208_366, 112_950, 321_316, 249_674, 71_642, 136_725]  # Sönderfallsenergierna för Lu-177.
Lu177_intensitet = [10.38, 6.2, 0.216, 0.2012, 0.1726, 0.047]  # Sönderfallsintensitet i % för respektive energi. Från laraweb.
Lu177_sannolikhet = np.cumsum(Lu177_intensitet) / np.sum(Lu177_intensitet)  # Kumulativa sannolikheten för sönderfall.

#   ----------------------------------------------------------------------
#   Filer med Data
#   ----------------------------------------------------------------------

tvärsnitt_file = '../given_data/Tvärsnittstabeller_Fotoner.xlsx'
attenueringsdata_file = '../given_data/Attenueringsdata.xlsx'
anatomidefinitioner_file = '../given_data/Anatomidefinitioner.xlsx'

mat_file = 'phantom_data.mat'


#   ----------------------------------------------------------------------
#   Funktioner
#   ----------------------------------------------------------------------

def plot_stuff(x_data, y_data, scatter, label_data,
               marker='o', color='blue', x_label='x-label', y_label='y-label', title='1',
               fig_size=(10, 10), symbol_size=100, font_size=30, alpha=1, line_width=10, x_lim=(0, 0), y_lim=(0, 0),
               grid=False, x_scale='linear', y_scale='linear'):
    """
    Funktion som skapar 2D plottar. Används nog inte i koden.
    """
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


def end_time(start):
    """
    Funktion som tar tiden för hur lång tid beräkningarna tar.
    """
    end = time.time()
    runtime = round((end - start), 1)
    if runtime < 60:
        print(f'Runtime: {runtime} seconds')
    else:
        print('Runtime: ' + str(round((runtime / 60), 1)) + ' minutes')


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