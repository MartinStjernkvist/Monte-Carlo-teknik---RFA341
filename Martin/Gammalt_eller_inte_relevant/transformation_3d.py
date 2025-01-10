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
import json

#   ----------------------------------------------------------------------
#   KONSTANTER
#   ----------------------------------------------------------------------

pi = np.pi
E_e = 0.511 * 10 ** 6  # eV
r_e = np.sqrt(0.07941)  # sqrt(b): re2 = e4/Ee2 ≈ 0.07941 b, https://en.wikipedia.org/wiki/Gamma_ray_cross_section
a_0 = 5.29177210903 * 10 ** (-11) * 10 ** (14)  # sqrt(b), bohr radius of hydrogen
c = 3 * 10 ** 8

radie_alpha = 1.2 * 10 ** (-15) * 4 ** (1 / 3)  # Radie alfapartikel i meter (Physics Handbook)

# enhet u
massa_H = 1
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
Lu177_intensitet = [10.38, 6.2, 0.216, 0.2012, 0.1726,
                    0.047]  # Sönderfallsintensitet i % för respektive energi. Från laraweb.
Lu177_sannolikhet = np.cumsum(Lu177_intensitet) / np.sum(Lu177_intensitet)  # Kumulativa sannolikheten för sönderfall.

At211_energi = [5.869, 5.2119, 5.1403, 4.9934, 4.895]  # Energi i MeV
At211_intensitet = [41.78, 0.0039, 0.0011, 0.0004, 0.00004]  # Intensitet i % för energierna
At211_sannolikhet = np.cumsum(At211_intensitet / np.sum(At211_intensitet))

#   ----------------------------------------------------------------------
#   Filer med Data
#   ----------------------------------------------------------------------

# tvärsnitt_file = '../given_data/Tvärsnittstabeller_Fotoner.xlsx'
# attenueringsdata_file = '../given_data/Attenueringsdata.xlsx'
# anatomidefinitioner_file = '../given_data/Anatomidefinitioner.xlsx'


tvärsnitt_file = 'Tvärsnittstabeller_Fotoner.xlsx'
attenueringsdata_file = 'Attenueringsdata.xlsx'
anatomidefinitioner_file = 'Anatomidefinitioner.xlsx'

mat_file = 'phantom_data.mat'


#Inläsning för vscode:
"""
tvärsnitt_file = r'given_data/Tvärsnittstabeller_Fotoner.xlsx'
attenueringsdata_file = r"given_data/Attenueringsdata.xlsx"
anatomidefinitioner_file = r"given_data/Anatomidefinitioner.xlsx"

mat_file = r"Martin/phantom_data.mat"
"""


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




def transformera_koordinatsystem(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B):
    """
    - Börjar på position A[x,y,z]- kalla detta koordinatsystem A.
    - Tar ett steg med steglängd steg_A_B, riktning (phi_A, theta_A), enligt koordinatsystemet i A.
            (exempelvis kommer phi_A att vara relativt enhetsvektorn i x-led för koord-syst A)
    - Tar ett steg till ny punkt - kalla denna punkt B.
    - Transformerar koordinatsystemet så att riktningsvektorn sammandfaller med
    nya koordinatsystemet.
            (nya enhetsvektorn i x-led, i B's koord-syst, ska ha samma riktning som fotonen
            hade när den tog steget)
    - Detta görs för att kunna sampla en ny riktning i nästa växelverkanprocess,
    då behövs nämligen ett koordinatsystem i B.

    :param steg_A_B: magnitud på steg från A till B
    :param phi_A: vinkel för steget mellan A och B
    :param theta_A: vinkel för steget mellan A och B

    :param steg_B_C: magnitud på steg från B till C
    :param phi_B: vinkel för steget mellan B och C
    :param theta_B: vinkel för steget mellan B och C

    :return: 1) vektor vars första 3 värden är positionen för punkt C enligt A's koord-syst
    """

    dx_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    dy_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    dz_A_B = steg_A_B * np.cos(theta_A)

    dx_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    dy_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    dz_B_C = steg_B_C * np.cos(theta_B)

    enhets_vektorer_A = np.array(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 0]
        ])

    # transformera med rotation i z-led (phi)
    # 4D matrix för att kunna skapa homogenitet-transformsmatris
    # R_z = np.array(
    #     [
    #         [np.cos(phi), -np.sin(phi), 0, 0],
    #         [np.sin(phi), np.cos(phi), 0, 0],
    #         [0, 0, 1, 0],
    #         [0, 0, 0, 0]
    #     ])

    R_z = np.array(
        [
            [np.cos(phi_A), -np.sin(phi_A), 0],
            [np.sin(phi_A), np.cos(phi_A), 0],
            [0, 0, 1]
        ])

    # R_y = np.array(
    #     [
    #         [np.cos(pi / 2 - theta), 0, -np.sin(pi / 2 - theta), 0],
    #         [0, 1, 0, 0],
    #         [np.sin(pi / 2 - theta), 0, np.cos(pi / 2 - theta), 0],
    #         [0, 0, 0, 0]
    #     ])

    # för att x-axeln ska sammanfalla med riktningsvektorn måste rotationsvinkeln vara pi/2 - theta
    angle = pi / 2 - theta_A
    R_y = np.array(
        [
            [np.cos(angle), 0, -np.sin(angle)],
            [0, 1, 0],
            [np.sin(angle), 0, np.cos(angle)],
        ])

    # R_y = np.identity(3)
    # R_z = np.identity(3)

    # R = R_y @ R_z

    # först rotation i theta (y-axeln), sedan rotation i phi (z-axeln)
    R = R_z @ R_y

    # print(f'\nR matrix: \n{R}')
    # print(f'\nlength of vector: \nx: {np.linalg.norm(R[0:3,0])}\ny: {np.linalg.norm(R[0:3,1])}\nz: {np.linalg.norm(R[0:3,2])}')

    Homogenous_matrix = np.array(
        [
            [R[0, 0], R[0, 1], R[0, 2], dx_A_B],
            [R[1, 0], R[1, 1], R[1, 2], dy_A_B],
            [R[2, 0], R[2, 1], R[2, 2], dz_A_B],
            [0, 0, 0, 1]
        ])

    vektor_A_C = Homogenous_matrix @ np.array(
        [
            [dx_B_C],
            [dy_B_C],
            [dz_B_C],
            [1]
        ])

    return vektor_A_C, dx_B_C, dy_B_C, dz_B_C


def transformera_koordinatsystem_extended(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B):
    dx_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    dy_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    dz_A_B = steg_A_B * np.cos(theta_A)

    dx_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    dy_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    dz_B_C = steg_B_C * np.cos(theta_B)

    enhets_vektorer_A = np.array(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 0]
        ])

    # transformera med rotation i z-led (phi)
    # 4D matrix för att kunna skapa homogenitet-transformsmatris
    # R_z = np.array(
    #     [
    #         [np.cos(phi), -np.sin(phi), 0, 0],
    #         [np.sin(phi), np.cos(phi), 0, 0],
    #         [0, 0, 1, 0],
    #         [0, 0, 0, 0]
    #     ])

    R_z = np.array(
        [
            [np.cos(phi_A), -np.sin(phi_A), 0],
            [np.sin(phi_A), np.cos(phi_A), 0],
            [0, 0, 1]
        ])

    # R_y = np.array(
    #     [
    #         [np.cos(pi / 2 - theta), 0, -np.sin(pi / 2 - theta), 0],
    #         [0, 1, 0, 0],
    #         [np.sin(pi / 2 - theta), 0, np.cos(pi / 2 - theta), 0],
    #         [0, 0, 0, 0]
    #     ])

    # för att x-axeln ska sammanfalla med riktningsvektorn måste rotationsvinkeln vara pi/2 - theta
    angle = pi / 2 - theta_A
    R_y = np.array(
        [
            [np.cos(angle), 0, -np.sin(angle)],
            [0, 1, 0],
            [np.sin(angle), 0, np.cos(angle)],
        ])

    # R_y = np.identity(3)
    # R_z = np.identity(3)

    # R = R_y @ R_z

    # först rotation i theta (y-axeln), sedan rotation i phi (z-axeln)
    R = R_z @ R_y

    # print(f'\nR matrix: \n{R}')
    # print(f'\nlength of vector: \nx: {np.linalg.norm(R[0:3,0])}\ny: {np.linalg.norm(R[0:3,1])}\nz: {np.linalg.norm(R[0:3,2])}')

    Homogenous_matrix = np.array(
        [
            [R[0, 0], R[0, 1], R[0, 2], dx_A_B],
            [R[1, 0], R[1, 1], R[1, 2], dy_A_B],
            [R[2, 0], R[2, 1], R[2, 2], dz_A_B],
            [0, 0, 0, 1]
        ])

    vektor_A_C = Homogenous_matrix @ np.array(
        [
            [dx_B_C],
            [dy_B_C],
            [dz_B_C],
            [1]
        ])

    enhets_vektorer_B = Homogenous_matrix @ enhets_vektorer_A

    return vektor_A_C, enhets_vektorer_B, dx_B_C, dy_B_C, dz_B_C


# @jit(nopython=True)
def ny_transformera_koordinatsystem(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B):
    """
    - Börjar på position A[x,y,z]- kalla detta koordinatsystem A.
    - Tar ett steg med steglängd steg_A_B, riktning (phi_A, theta_A), enligt koordinatsystemet i A.
            (exempelvis kommer phi_A att vara relativt enhetsvektorn i x-led för koord-syst A)
    - Tar ett steg till ny punkt - kalla denna punkt B.
    - Transformerar koordinatsystemet så att riktningsvektorn sammandfaller med
    nya koordinatsystemet.
            (nya enhetsvektorn i x-led, i B's koord-syst, ska ha samma riktning som fotonen
            hade när den tog steget)
    - Detta görs för att kunna sampla en ny riktning i nästa växelverkanprocess,
    då behövs nämligen ett koordinatsystem i B.

    :param steg_A_B: magnitud på steg från A till B
    :param phi_A: vinkel för steget mellan A och B
    :param theta_A: vinkel för steget mellan A och B

    :param steg_B_C: magnitud på steg från B till C
    :param phi_B: vinkel för steget mellan B och C
    :param theta_B: vinkel för steget mellan B och C

    :return: 3 värden är förflyttningen från B till C, enligt A's koord-syst
    """

    dx_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    dy_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    dz_A_B = steg_A_B * np.cos(theta_A)

    dx_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    dy_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    dz_B_C = steg_B_C * np.cos(theta_B)

    # transformera med rotation i z-led (phi)
    R_z = np.array(
        [
            [np.cos(phi_A), -np.sin(phi_A), 0],
            [np.sin(phi_A), np.cos(phi_A), 0],
            [0, 0, 1]
        ], dtype=np.float64)

    # för att x-axeln ska sammanfalla med riktningsvektorn måste rotationsvinkeln vara pi/2 - theta
    angle = pi / 2 - theta_A
    R_y = np.array(
        [
            [np.cos(angle), 0, -np.sin(angle)],
            [0, 1, 0],
            [np.sin(angle), 0, np.cos(angle)],
        ], dtype=np.float64)

    # först rotation i theta (y-axeln), sedan rotation i phi (z-axeln)
    R = R_z @ R_y

    Homogenous_matrix = np.array(
        [
            [R[0, 0], R[0, 1], R[0, 2], dx_A_B],
            [R[1, 0], R[1, 1], R[1, 2], dy_A_B],
            [R[2, 0], R[2, 1], R[2, 2], dz_A_B],
            [0, 0, 0, 1]
        ], dtype=np.float64)

    vektor_A_C = Homogenous_matrix @ np.array(
        [
            [dx_B_C],
            [dy_B_C],
            [dz_B_C],
            [1]   # nödvändigt för beräkningen
        ], dtype=np.float64)

    vektor_A_B = np.array([
        [dx_A_B],
        [dy_A_B],
        [dz_A_B],
        [1]  # nödvändigt för beräkningen
    ], dtype=np.float64)

    # Vill ha vektor B->C
    vektor = np.subtract(vektor_A_C, vektor_A_B)
    dx, dy, dz = vektor[0][0] / voxel_sidlängd, vektor[2][0] / voxel_sidlängd, vektor[3][0] / voxel_sidlängd

    return dx, dy, dz


if __name__ == "__main__":
    #   ----------------------------------------------------------------------
    #   INPUT
    #   ----------------------------------------------------------------------

    # startpunkt A
    x_start = 1
    y_start = 1
    z_start = 1

    x_A = x_start
    y_A = y_start
    z_A = z_start

    # steg 1: från A till B
    theta_A = 3 * pi / 8
    phi_A = 1 * pi / 8
    steg_A_B = 3

    # steg 2: Från B till C, enligt koordinatsystemet för B
    theta_B = pi / 3
    phi_B = pi / 3
    steg_B_C = 2

    #   ----------------------------------------------------------------------
    #   BERÄKNING
    #   ---------------------------------------------------------------------

    dx_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    dy_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    dz_A_B = steg_A_B * np.cos(theta_A)

    vektor_A_C, enhets_vektorer_B, dx_B_C, dy_B_C, dz_B_C = transformera_koordinatsystem_extended(steg_A_B, phi_A,
                                                                                                  theta_A,
                                                                                                  steg_B_C, phi_B,
                                                                                                  theta_B)

    x_A_C = vektor_A_C[0]
    x_tot = x_A + x_A_C

    y_A_C = vektor_A_C[1]
    y_tot = y_A + y_A_C

    z_A_C = vektor_A_C[2]
    z_tot = z_A + z_A_C

    #   ----------------------------------------------------------------------
    #   FELSÖKNING
    #   ---------------------------------------------------------------------

    print(f'tot: {np.array([[x_tot], [y_tot], [z_tot]])}')

    print(f'\nresult: \n{vektor_A_C}\nenhets_vektorer_B: \n{enhets_vektorer_B}')

    print(
        f'\ndot product x ({enhets_vektorer_B[0:3, 0]}), y ({enhets_vektorer_B[0:3, 1]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 1])}')
    print(
        f'dot product y ({enhets_vektorer_B[0:3, 1]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 1], enhets_vektorer_B[0:3, 2])}')
    print(
        f'dot product x ({enhets_vektorer_B[0:3, 0]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 2])}')

    #   ----------------------------------------------------------------------
    #   PLOTTNING
    #   ---------------------------------------------------------------------

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection='3d')
    ax.set_xlim([0, 6])
    ax.set_ylim([0, 6])
    ax.set_zlim([0, 6])

    ax.scatter(x_start, y_A, z_A, label='A', color='black', s=100)
    ax.scatter(x_start + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, label='B', color='grey', s=100)
    ax.scatter(x_A + x_A_C, y_A + y_A_C, z_A + z_A_C, label='C', color='purple', s=100)

    ax.quiver(x_A, y_A, z_A, 1, 0, 0, color='orange')
    ax.quiver(x_A, y_A, z_A, 0, 1, 0, color='brown')
    ax.quiver(x_A, y_A, z_A, 0, 0, 1, color='green')

    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 0], enhets_vektorer_B[1, 0],
              enhets_vektorer_B[2, 0],
              color='orange')
    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 1], enhets_vektorer_B[1, 1],
              enhets_vektorer_B[2, 1],
              color='brown')
    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 2], enhets_vektorer_B[1, 2],
              enhets_vektorer_B[2, 2],
              color='green')

    ax.quiver(x_A, y_A, z_A, dx_A_B, dy_A_B, dz_A_B, color='blue', label='första steget, A-> B')
    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, (x_A_C - dx_A_B), (y_A_C - dy_A_B), (z_A_C - dz_A_B),
              color='magenta',
              label='andra steget, B->C')
    ax.quiver(x_A, y_A, z_A, x_A_C, y_A_C, z_A_C, color='red')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.legend()
    plt.show()

    #   ----------------------------------------------------------------------
    #   INPUT - iteration 2
    #   ---------------------------------------------------------------------

    # startpunkt
    x_B = x_A + dx_A_B
    y_B = y_A + dy_A_B
    z_B = z_A + dz_A_B

    # steg 2
    # steg_B_C
    # phi_B
    # theta_B

    # steg 3
    steg_C_D = 2
    theta_C = pi / 4
    phi_C = pi / 4

    #   ----------------------------------------------------------------------
    #   BERÄKNING
    #   ---------------------------------------------------------------------

    dx_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    dy_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    dz_B_C = steg_B_C * np.cos(theta_B)

    vektor_B_D, enhets_vektorer_C, dx_C_D, dy_C_D, dz_C_D = transformera_koordinatsystem_extended(steg_B_C, phi_B,
                                                                                                  theta_B,
                                                                                                  steg_C_D, phi_C,
                                                                                                  theta_C)

    x_B_D = vektor_B_D[0]
    x_tot = x_B + x_B_D

    y_B_D = vektor_B_D[1]
    y_tot = y_B + y_B_D

    z_B_D = vektor_B_D[2]
    z_tot = z_B + z_B_D

    #   ----------------------------------------------------------------------
    #   FELSÖKNING
    #   ---------------------------------------------------------------------

    # print(f'tot: {np.array([[x_tot], [y_tot], [z_tot]])}')
    #
    # print(f'\nresult: \n{vektor_A_C}\nenhets_vektorer_B: \n{enhets_vektorer_B}')
    #
    # print(
    #     f'\ndot product x ({enhets_vektorer_B[0:3, 0]}), y ({enhets_vektorer_B[0:3, 1]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 1])}')
    # print(
    #     f'dot product y ({enhets_vektorer_B[0:3, 1]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 1], enhets_vektorer_B[0:3, 2])}')
    # print(
    #     f'dot product x ({enhets_vektorer_B[0:3, 0]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 2])}')

    #   ----------------------------------------------------------------------
    #   PLOTTNING
    #   ---------------------------------------------------------------------

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection='3d')
    ax.set_xlim([0, 6])
    ax.set_ylim([0, 6])
    ax.set_zlim([0, 6])

    ax.scatter(x_B, y_B, z_B, label='B', color='black', s=100)
    ax.scatter(x_B + dx_B_C, y_B + dy_B_C, z_B + dz_B_C, label='C', color='grey', s=100)
    ax.scatter(x_B + x_B_D, y_B + y_B_D, z_B + z_B_D, label='D', color='purple', s=100)

    ax.quiver(x_B + dx_B_C, y_B + dy_B_C, z_B + dz_B_C, enhets_vektorer_C[0, 0], enhets_vektorer_C[1, 0],
              enhets_vektorer_C[2, 0],
              color='orange')
    ax.quiver(x_B + dx_B_C, y_B + dy_B_C, z_B + dz_B_C, enhets_vektorer_C[0, 1], enhets_vektorer_C[1, 1],
              enhets_vektorer_C[2, 1],
              color='brown')
    ax.quiver(x_B + dx_B_C, y_B + dy_B_C, z_B + dz_B_C, enhets_vektorer_C[0, 2], enhets_vektorer_C[1, 2],
              enhets_vektorer_C[2, 2],
              color='green')

    ax.quiver(x_B, y_B, z_B, dx_B_C, dy_B_C, dz_B_C, color='blue', label='första steget, B -> C')
    ax.quiver(x_B + dx_B_C, y_B + dy_B_C, z_B + dz_B_C, (x_B_D - dx_B_C), (y_B_D - dy_B_C), (z_B_D - dz_B_C),
              color='magenta',
              label='andra steget, C -> D')
    ax.quiver(x_B, y_B, z_B, x_B_D, y_B_D, z_B_D, color='red')
    ax.quiver(0, 0, 0, x_tot, y_tot, z_tot, color='black')

    ax.scatter(x_start, y_A, z_A, label='A', color='black', s=100)

    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 0], enhets_vektorer_B[1, 0],
              enhets_vektorer_B[2, 0],
              color='orange')
    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 1], enhets_vektorer_B[1, 1],
              enhets_vektorer_B[2, 1],
              color='brown')
    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 2], enhets_vektorer_B[1, 2],
              enhets_vektorer_B[2, 2],
              color='green')

    ax.quiver(x_A, y_A, z_A, dx_A_B, dy_A_B, dz_A_B, color='blue', label='första steget, A-> B')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.legend()
    plt.show()
