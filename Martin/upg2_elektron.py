"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IMPORTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""

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

r_e_m = 2.81 * 10 ** (-15)  # Elektronradie i meter.

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

At211_energi_MeV = [5.869, 5.2119, 5.1403, 4.9934, 4.895]  # Sönderfallsdata för At-211.
At211_energi = [(lambda x: x * 10 ** 6)(x) for x in At211_energi_MeV]  # Energi i eV
At211_intensitet = [41.78, 0.0039, 0.0011, 0.0004, 0.00004]  # Intensitet i % för energierna
At211_sannolikhet = np.cumsum(At211_intensitet) / np.sum(At211_intensitet)

rho_vatten = 998  # Vattens densitet i kg / m^3.

#   ----------------------------------------------------------------------
#   Filer med Data
#   ----------------------------------------------------------------------

# Uppgift 2: Elektroner.
Y90_file = '../given_data/Y90_Spektrum.xlsx'
elektron_stopping_power_data = 'Elektron_stopping_power_range_data'
elektron_scatter_power_data = 'Elektron_scatter_power_vatten_data'

# Inläsning för vscode:
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STARTENERGI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""

#   ----------------------------------------------------------------------
#   Sönderfallsdata för Y-90.
#   ----------------------------------------------------------------------

file_Y90 = pd.read_excel(Y90_file)

# Ta ut mätpunkter med intensitet och energi.
Energi_Y90 = file_Y90['Energy (MeV)']  # MeV
Intensitet_Y90 = file_Y90['#/nt']


#   ----------------------------------------------------------------------
#   Kurvanpassning till sönderfallsdatan.
#   ----------------------------------------------------------------------

def polynom_funktion(x, a, b, c, d, e, f):
    return a * x ** 5 + b * x ** 4 + c * x ** 3 + d * x ** 2 + e * x + f


params, cv = curve_fit(polynom_funktion, Energi_Y90, Intensitet_Y90)
a, b, c, d, e, f = params

olika_energier = np.linspace(np.min(Energi_Y90), np.max(Energi_Y90), 10_000)


#   ----------------------------------------------------------------------
#   Använd rejektionsmetoden för att sampla elektronenergi.
#   ----------------------------------------------------------------------
def elektron_energi_start():
    elektron_energi = 0
    hittat = 0
    while hittat == 0:
        x_rand = np.random.rand() * np.max(Energi_Y90)
        y_rand = np.random.rand() * np.max(Intensitet_Y90)

        if y_rand < polynom_funktion(x_rand, *params):
            hittat = 1
            elektron_energi = x_rand

        else:
            hittat = 0

    return elektron_energi


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STARTRIKTNING

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


@jit(nopython=True)
def riktning_uniform():
    """
    Uniform riktning.
    """
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi


@jit(nopython=True)
def riktning_skal():
    """
    Ifall skalförfördelnning:
    pi / 2 < phi < 3 * pi / 2 för att effektivisera koden.
    Då kan antalet iterationer halveras, eftersom man vet att
    hälften av partiklarna ändå skulle lämna sfären.
    """
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = pi / 2 * (2 * np.random.rand() + 1)
    return theta, phi


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STARTPOSITION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


@jit(nopython=True)
def position_start_innanför(radie_sfär):
    """
    Funktion som samplar startposition.
    Fördelningsfunktion: uniform innanför sfär.
    """
    r = radie_sfär * np.random.rand()

    x = r
    # Utnyttja sfärisk symmetri.
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
def position_start_skal(radie_sfär, radie_partikel):
    """
    Funktion som samplar startposition.
    Fördelningsfunktion: ytfördelning på sfär.
    """
    r = radie_sfär - 0.5 * radie_partikel  # För att inte endast theta = pi ska ge utslag

    x = r
    # Utnyttja sfärisk symmetri.
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELEKTRON POLARVINKEL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def polar_vinkel(steglängd, scattering_power):
    """
    Funktion som ger polarvinkeln theta för elektronspridning.
    :param steglängd: Steglängden innan spridningen.
    :param scattering_power: Scattering power från annan funktion.
    :return: Spridningsvinkeln theta (sfäriska koordinater).
    """

    # steglängd_cm = steglängd * 10 ** 2

    R = np.random.random()  # Slumpmässig tal mellan 0-1

    T = scattering_power

    theta_ny = np.sqrt(-T * steglängd * np.log(1 - R))

    print(f'polarvinkel: {theta_ny * 360 / (2 * np.pi):.3f} grader')
    return theta_ny


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCATTERING POWER FRÅN ENERGI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def scattering_power_från_energi(elektron_energi, scatterpower_data, rho_medium):
    """
    Funktion som linjärinterpolerar scattering power utifrån energi.
    :param elektron_energi: Elektronens energi.
    :param scatterpower_data: Scattering power tabelldata.
    :param rho_medium: Mediumets densitet (kg/m^3).
    :return: Scattering power (rad^2 / cm).
    """

    # Omvandlar från MeV till eV och till en np.array
    energi_MeV_list = scatterpower_data[:, 0]
    energi_list = np.array([(lambda x: x * 10 ** 6)(x) for x in energi_MeV_list])

    mass_scattering_power_list = np.array(scatterpower_data[:, 1])

    # Hittar närmast energi som liknar elektron energin och ta fram scatter power:

    # Tar index för närmaste energi på elektronen.
    diff = np.abs(energi_list - elektron_energi)
    closest_indices = np.argsort(diff)[:2]

    # Få scattering power nära energierna och energivärdena närmast.
    T_close_cm = mass_scattering_power_list[closest_indices] * rho_medium * 10 ** (-3)  # Densitet: kg/m^3 till g/cm^3
    T_close = T_close_cm * 10 ** 2  # scattering power / m

    energi_close = energi_list[closest_indices]

    # Linjärinterpolering av scattering power.
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        T = T_close[0]

    else:
        T = T_close[0] + (elektron_energi - energi_close[0]) * (
                T_close[1] - T_close[0]) / (
                    energi_close[1] - energi_close[0])

    print(f'scattering power T / cm: {T / 100}')

    return T


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STOPPING POWER OCH STEGLÄNGD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def stopping_power_och_steglängd(energi, rho_medium, stopping_power_data):
    """
    Funktion som ger stopping power och
    :param energi: Partikelns energi i eV.
    :param rho_medium: Mediumets densitet i kg / m^3.
    :param stopping_power_data: Tabellerad data.
    :return: Stopping power (eV/m) och steglängden (m).
    """

    energi_MeV_list = stopping_power_data[:, 0]  # MeV
    energi_list = np.array([(lambda x: x * 10 ** 6)(x) for x in energi_MeV_list])  # eV
    # print(f'energi list: {energi_list}')

    stopping_power_MeV_cm_g_list = stopping_power_data[:, 1]  # MeV cm^2 / g
    stopping_power_list = np.array(
        [(lambda x: x * 10 ** (6 - 2 * 2 + 3))(x) for x in stopping_power_MeV_cm_g_list])  # eV m^2 / kg
    # print(f'stopping power list: {stopping_power_list}')

    CSDA_cm_g_list = stopping_power_data[:, 2]  # g / cm^2
    CSDA_list = np.array([(lambda x: x * 10 ** (-3 + 2 * 2))(x) for x in CSDA_cm_g_list])  # kg / m^2
    # print(f'CSDA list: {CSDA_list}')

    # Tar index för närmaste energi på alfapartikeln.
    diff = np.abs(energi_list - energi)
    closest_indices = np.argsort(diff)[:2]

    # Närmsta värdena:
    energi_close = energi_list[closest_indices]
    stopping_power_close = stopping_power_list[closest_indices]
    CSDA_close = CSDA_list[closest_indices]

    # Linjärinterpolera och få fram stopping power och slumpmässig steglängd.
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        stopping_power = stopping_power_close[0]
        CSDA = CSDA_close[0]

    else:
        stopping_power = (stopping_power_close[0] + (energi - energi_close[0]) * (
                stopping_power_close[1] - stopping_power_close[0]) / (
                                  energi_close[1] - energi_close[0]))

        CSDA = (CSDA_close[0] + (energi - energi_close[0]) * (
                CSDA_close[1] - CSDA_close[0]) / (
                        energi_close[1] - energi_close[0]))

    steglängd = CSDA / rho_medium
    stopping_power = stopping_power * rho_medium  # eV / m

    return stopping_power, steglängd


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STEGLÄNGD FRÅN ENERGI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def steglängd_från_energi(energi_start, rho_medium, stopping_power_data, energiförlust_faktor):
    """
    Funktion som beräknar steglängden utifrån startenergi och energiförlust.
    :param energi_start: Startenergin för elektronen.
    :param rho_medium: Mediumets densitet (kg/m^3).
    :param stopping_power_data: Stopping power tabelldata.
    :param energiförlust_faktor: Faktor för energiförlust, t.ex. 0.97 för 3%.
    :return: Steglängden till nästa kollision.
    """

    energi_1 = energi_start

    # Bestäm ny energi, varefter steglängden kan lösas ut m.h.a. stopping power.
    energi_2 = energi_1 * energiförlust_faktor  # antag energiförlust innan nästa kollission

    # Stopping power från tabelldata, med hjälp av linjärinterpolering.
    # Delvis stopping power för energi_1, delvis för energi_2.
    stopping_power_1, _ = stopping_power_och_steglängd(energi_1, rho_medium, stopping_power_data)
    stopping_power_2, _ = stopping_power_och_steglängd(energi_2, rho_medium, stopping_power_data)

    # Medelvärdet av stopping power mellan startpunkt och slutpunkt.
    stopping_power_medel = (stopping_power_2 + stopping_power_1) / 2
    print(f'stopping power medel: {stopping_power_medel * 10 ** (-6):.3f} MeV / m')

    # Beräkna steglängden.
    steglängd = - (energi_2 - energi_1) / stopping_power_medel
    print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

    return steglängd


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TRANSFORMERA KOORDINATSYSTEM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


# @jit(nopython=True)
def ny_steg_transformera_koordinatsystem_3d(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B):
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

    # för att z-axeln ska sammanfalla med riktningsvektorn måste rotationsvinkeln vara theta_A
    angle = theta_A
    R_y = np.array(
        [
            [np.cos(angle), 0, np.sin(angle)],
            [0, 1, 0],
            [-np.sin(angle), 0, np.cos(angle)],
        ], dtype=np.float64)

    # först rotation i theta (y-axeln), sedan rotation i phi (z-axeln)
    R = np.dot(R_z, R_y)

    Homogenous_matrix = np.eye(4, dtype=np.float64)
    Homogenous_matrix[:3, :3] = R
    Homogenous_matrix[:3, 3] = np.array([dx_A_B, dy_A_B, dz_A_B], dtype=np.float64)

    vektor_A_C = np.dot(Homogenous_matrix, np.array([dx_B_C, dy_B_C, dz_B_C, 1.0], dtype=np.float64))

    vektor_A_B = np.array([dx_A_B, dy_A_B, dz_A_B, 1.0], dtype=np.float64)

    # Vill ha vektor B->C
    vektor = vektor_A_C - vektor_A_B
    dx, dy, dz = vektor[0], vektor[1], vektor[2]

    return dx, dy, dz


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LADDAD PARTIKEL VÄG - ELEKTRONER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def laddad_partikel_väg_elektron(energi_start, position_start, phi, theta, radie_sfär, rho_medium,
                                 stopping_power_data, scatter_power_data, energiförlust_faktor):
    """
    Funktion som följer alfapartikeln allteftersom den växelverkar i ett medium.
    :param energi_start: Elektronens startenergi.
    :param position_start: Elektronens startposition.
    :param phi: Sfärisk vinkel.
    :param theta: Sfärisk vinkel.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param rho_medium: Mediumets densitet (kg / m^3).
    :param stopping_power_data: Stopping power tabelldata.
    :param scatter_power_data: Scatter power tabelldata.
    :param energiförlust_faktor: Energiförlustfaktor efter varje steglängd.
    :return: Energideponeringen innanför sfären.
    """

    # Initiera tomma listor för att spara datan.
    x, y, z, dos = [], [], [], []

    # Startposition och startenergi.
    position_vektor = position_start
    energi = energi_start

    # Steglängd tas så att en viss energiförlust sker.
    steglängd = steglängd_från_energi(energi, rho_medium, stopping_power_data, energiförlust_faktor)

    # Den nya energin för elektronen.
    energi_ny = energi * energiförlust_faktor

    # Första steget.
    dx = steglängd * np.sin(theta) * np.cos(phi)
    dy = steglängd * np.sin(theta) * np.sin(phi)
    dz = steglängd * np.cos(theta)

    position_vektor += np.array([dx, dy, dz])
    # print(f'positionsvektor: {position_vektor}')

    if np.sqrt(np.dot(position_vektor, position_vektor)) < radie_sfär:
        dos.append(energi - energi_ny)
        x.append(position_vektor[0])
        y.append(position_vektor[1])
        z.append(position_vektor[2])

    energi = energi_ny

    #   ----------------------------------------------------------------------
    #   Medan elektronen befinner sig i sfären
    #   och
    #   har en energi som är över en tröskelenergi.
    #   ----------------------------------------------------------------------
    while np.sqrt(np.dot(position_vektor, position_vektor)) < radie_sfär and energi > energi_start * 10 ** (-6):

        # Steglängd tas så att en viss energiförlust sker.
        steglängd_ny = steglängd_från_energi(energi, rho_medium, stopping_power_data, energiförlust_faktor)

        # Den nya energin för elektronen.
        energi_ny = energi * energiförlust_faktor

        # Scattering power för den nya energin, efter steget.
        scattering_power = scattering_power_från_energi(energi_ny, scatter_power_data, rho_medium)

        # Polarvinkel utifrån scattering power.
        theta_ny = polar_vinkel(steglängd, scattering_power)
        phi_ny = np.random.random() * 2 * pi

        # Ändrar på positionsvektor efter att transformations matrisen
        dx, dy, dz = ny_steg_transformera_koordinatsystem_3d(steglängd, phi, theta, steglängd_ny, phi_ny, theta_ny)
        position_vektor += np.array([dx, dy, dz])

        #   ----------------------------------------------------------------------
        #   Håll reda på ifall partikeln befinner sig i sfären eller inte.
        #   ----------------------------------------------------------------------
        if np.sqrt(np.dot(position_vektor, position_vektor)) > radie_sfär:
            break

        # Om elektronen fortfarande är i sfären -> spara mätpunkter.
        else:
            # Mätpunkter för att plotta dosfördelningen.
            dos.append(energi - energi_ny)
            x.append(position_vektor[0])
            y.append(position_vektor[1])
            z.append(position_vektor[2])

            # Ändrar på vinklarna, värdet på steglängden och energin innan nästa vinkeländring.
            phi = phi_ny
            theta = theta_ny
            steglängd = steglängd_ny
            energi = energi_ny

            print(f'energi: {energi:.3f} eV')
            continue

    # Beräkna den totala energideponeringen (i sfären) från partikeln.
    energideponering = energi_start - energi
    # print(f'energideponering: {energideponering} eV')
    return energideponering, x, y, z, dos


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIMULERINGSKOD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                    position_start, radie_sfär, energiförlust_faktor):
    """
    Monte-Carlo simulering för elektrontransport.
    :param iterationer: Antal sönderfall som ska simuleras.
    :param rho_medium: Mediumets densitet (kg / m^3).
    :param radie_partikel: Partikelns radie.
    :param stopping_power_data: Stopping power tabelldata.
    :param scatter_power_data: Scatter power tabelldata.
    :param position_start: Uniform fördelning i sfären, eller ytfördelning.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param energiförlust_faktor: Energiförlustfaktor efter varje steglängd.
    :return: Summeringen av energideponeringen innanför sfären.
    """

    # Initiera en energisumma och tomma listor för att spara datan.
    energideponering_summa = 0
    x_list, y_list, z_list, dos_list = [], [], [], []

    #   ----------------------------------------------------------------------
    #   Vilken fördelningsfunktion som ska användas bestämmer hur
    #   sampling av riktning och position sker.
    #   ----------------------------------------------------------------------
    if position_start == position_start_skal:
        #   ----------------------------------------------------------------------
        #   Ytfördelning på en sfär.
        #   ----------------------------------------------------------------------
        iterationer = 0.5 * iterationer
        for i in range(int(iterationer)):
            # Sampla riktning och startposition.
            theta, phi = riktning_skal()
            position_start = position_start_skal(radie_sfär, radie_partikel)

            # Sampla startenergin.
            energi_start = elektron_energi_start() * 10 ** 6  # i eV
            # print(f'energi: {energi_start * 10 ** (-6)} MeV')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering, x, y, z, dos = laddad_partikel_väg_elektron(energi_start, position_start, phi, theta,
                                                                          radie_sfär, rho_medium,
                                                                          stopping_power_data, scatter_power_data,
                                                                          energiförlust_faktor)

            # Spara mätpunkter för plottning.
            x_list += x
            y_list += y
            z_list += z
            dos_list += dos

            # Summera energideponeringsbidragen från respektive iteration.
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    else:
        #   ----------------------------------------------------------------------
        #   Uniform fördelning i en sfär.
        #   ----------------------------------------------------------------------
        for i in range(iterationer):
            # Sampla startenergin.
            energi_start = elektron_energi_start() * 10 ** 6  # i eV
            # print(f'energi: {energi_start * 10 ** (-6)} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_uniform()
            position_start = position_start_innanför(radie_sfär)

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering, x, y, z, dos = laddad_partikel_väg_elektron(energi_start, position_start, phi, theta,
                                                                          radie_sfär, rho_medium,
                                                                          stopping_power_data, scatter_power_data,
                                                                          energiförlust_faktor)
            # Spara mätpunkter för plottning.
            x_list += x
            y_list += y
            z_list += z
            dos_list += dos

            # Summera energideponeringsbidragen från respektive iteration.
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    # print('total energideponering: ', energideponering_summa)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x_list, y_list, z_list, c=dos_list, cmap='plasma', label='Partikel position')
    # Fixa colorbar för att se energideponeringen i figuren

    # fig.colorbar(ax=ax, label='Energideponering',)

    ax.set_xlabel('x-axel (m)')
    ax.set_ylabel('y-axel (m)')
    ax.set_ylabel('z-axel (m)')

    # Testar att sätta en sfär för tumören
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = radie_sfär * np.cos(u) * np.sin(v)
    y = radie_sfär * np.sin(u) * np.sin(v)
    z = radie_sfär * np.cos(v)
    ax.plot_wireframe(x, y, z, color="k", alpha=0.3, label='Tumören')
    ax.legend()

    # Visa figur
    plt.show()

    return energideponering_summa


def energideponering_eV_till_Gy(energideponering_eV, rho_medium, radie_sfär):
    """
    Beräkna absorberad i termer av Gy.
    :param energideponering_eV: Absorberad dos i eV.
    :param rho_medium: Mediumets densitet.
    :param radie_sfär: Sfärens radie.
    """

    # Beräkna volym (m^3) och massa (kg).
    V = 4 / 3 * np.pi * radie_sfär ** 3
    massa = V * rho_medium

    # Beräkna absorberad dos i Joule.
    energideponering_J = energideponering_eV * 1.602 * 10 ** (-19)

    # Beräkna absorberad dos i Gy.
    energideponering_Gy = energideponering_J / massa

    return energideponering_Gy


if __name__ == "__main__":
    #   ----------------------------------------------------------------------
    #   Köra simuleringskoden.
    #   ----------------------------------------------------------------------
    iterationer = 10 ** 3
    dummy_iterationer = 10 ** 2

    energiförlust_faktor = 0.95  # Energi efter en kollision = 98% av föregående energin.

    stopping_power_data = np.loadtxt(elektron_stopping_power_data)
    scatter_power_data = np.loadtxt(elektron_scatter_power_data)

    rho_medium = rho_vatten
    radie_partikel = r_e_m

    radie_sfär_skal = 300 * 10 ** (-6)  # m
    radie_sfär_innanför = 1 * 10 ** (-3)  # m

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_elektron(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                        position_start_skal,
                        radie_sfär_skal, energiförlust_faktor)

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                scatter_power_data,
                                                position_start_skal,
                                                radie_sfär_skal, energiförlust_faktor)

    energideponering_skal_Gy = energideponering_eV_till_Gy(energideponering_tot_skal, rho_medium, radie_sfär_skal)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_elektron(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                        position_start_innanför,
                        radie_sfär_innanför, energiförlust_faktor)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                    scatter_power_data,
                                                    position_start_innanför, radie_sfär_innanför, energiförlust_faktor)

    energideponering_innanför_Gy = energideponering_eV_till_Gy(energideponering_tot_innanför, rho_medium,
                                                               radie_sfär_innanför)
    end_time(start)

    print(
        '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nRESULTAT\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')

    print(
        f'\nSkal (300 mikrometer): Energideponering:\n{energideponering_skal_Gy * 10 ** 8 / iterationer} E-08 Gy / sönderfall')
    print(f'faktor {(energideponering_skal_Gy * 10 ** 8 / iterationer) / 4.07} av facit')

    print(
        f'\nInnanför (1 mm): Energideponering:\n{energideponering_innanför_Gy * 10 ** 9 / iterationer} E-09 Gy / sönderfall')
    print(f'faktor {(energideponering_innanför_Gy * 10 ** 9 / iterationer) / 5.22} av facit')
