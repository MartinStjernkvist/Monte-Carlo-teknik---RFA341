"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IMPORTS
Pythonpaket och konstanter.
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

#   -----------------------------------
#   KONSTANTER
#   -----------------------------------

font_size = 20
font_size_title = font_size * 1.25

pi = np.pi
E_e = 0.511 * 10 ** 6  # eV
r_e = np.sqrt(0.07941)  # sqrt(b): re2 = e4/Ee2 ≈ 0.07941 b, https://en.wikipedia.org/wiki/Gamma_ray_cross_section
a_0 = 5.29177210903 * 10 ** (-11) * 10 ** 14  # sqrt(b), bohr radius of hydrogen
c = 3 * 10 ** 8

r_e_m = 2.81 * 10 ** (-15)  # Elektronradie i meter.

radie_alpha = 1.2 * 10 ** (-15) * 4 ** (1 / 3)  # Radie alfapartikel i meter (Physics Handbook)

At211_energi_MeV = [5.869, 5.2119, 5.1403, 4.9934, 4.895]  # Sönderfallsdata för At-211.
At211_energi = [(lambda x: x * 10 ** 6)(x) for x in At211_energi_MeV]  # Energi i eV
At211_intensitet = [41.78, 0.0039, 0.0011, 0.0004, 0.00004]  # Intensitet i % för energierna
At211_sannolikhet = np.cumsum(At211_intensitet) / np.sum(At211_intensitet)

rho_vatten = 998  # Vattens densitet i kg / m^3.

#   -----------------------------------
#   Filer med Data
#   -----------------------------------

# Uppgift 2: Alfapartiklar
stopping_power_alfa_file = '../Stoppingpower_data_alfa'

# Uppgift 2: Elektroner.
Y90_file = '../../given_data/Y90_Spektrum.xlsx'
elektron_stopping_power_data = 'Elektron_stopping_power_range_data'
elektron_scatter_power_data = 'Elektron_scatter_power_vatten_data'


#   -----------------------------------
#   Funktioner
#   -----------------------------------

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
STARTRIKTNING
Samplar startriktningen för en partikel, utifrån fördelningsfunktion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


@jit(nopython=True)
def riktning_uniform():
    """
    Uniform riktning i en sfär.
    """
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi


@jit(nopython=True)
def riktning_skal():
    """
    Ifall ytfördelning (skalet på en sfär):
    pi / 2 < phi < 3 * pi / 2 för att effektivisera koden.
    Då kan antalet iterationer halveras, eftersom man vet att
    hälften av partiklarna ändå skulle lämna sfären.
    """
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = pi / 2 * (2 * np.random.rand() + 1)
    return theta, phi

"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STOPPING POWER OCH STEGLÄNGD
Beräkna stopping power och steglängden utifrån energi.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def stopping_power_och_steglängd(energi, rho_medium, stopping_power_data):
    """
    Funktion som ger stopping power och
    :param energi: Partikelns energi (eV).
    :param rho_medium: Mediumets densitet (kg/m^3).
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
    stopping_power = stopping_power * rho_medium

    return stopping_power, steglängd


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STARTPOSITION
Samplar startposition för partiklarna utifrån fördelningsfunktion.
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
STARTENERGI
Beräknar startenergin för partikeln, utifrån söndefallsdata.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def energi_start(radionuklid_energi, radionuklid_sannolikhet):
    """
    Funktion som samplar den ursprungliga fotonenergin för varje ny foton som sänds ut.
    :param radionuklid_energi: Lista med fotonenergier som är möjliga.
    :param radionuklid_sannolikhet: Lista med sannolikheten för fotonergierna.
    :return:
    """

    slump_tal = np.random.rand()

    if slump_tal <= radionuklid_sannolikhet[0]:
        foton_energi = radionuklid_energi[0]
    elif slump_tal <= radionuklid_sannolikhet[1]:
        foton_energi = radionuklid_energi[1]
    elif slump_tal <= radionuklid_sannolikhet[2]:
        foton_energi = radionuklid_energi[2]
    elif slump_tal <= radionuklid_sannolikhet[3]:
        foton_energi = radionuklid_energi[3]
    elif slump_tal <= radionuklid_sannolikhet[4]:
        foton_energi = radionuklid_energi[4]
    else:
        foton_energi = radionuklid_energi[5]

    return foton_energi


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENERGIFÖRLUST
Beräknar energiförlusten efter varje steg som partikeln tar.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def energi_efter_energiförlust(energi, steglängd, rho_medium, stopping_power_data):
    """
    Funktion som beräknar energiförlusten efter varje steg som partikeln tar.
    :param energi: Partikelns energi innan steget.
    :param steglängd: Stegstorleken (m).
    :param rho_medium: Densiteten av mediumet (kg/m^3).
    :param stopping_power_data: Stopping power data.
    :return: Energiförlusten under steget.
    """

    # Beräkna stopping power, utifrån tabellerad data.
    stopping_power, _ = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
    # print(f'stopping power: {stopping_power}')

    # Beräkna energiförlusten med hjälp av stopping power.
    energiförlust = stopping_power * steglängd
    # print(f'energiförlust: {energiförlust * 10**(-6):.4f} MeV')

    # Beräkna partikelns nya energi, efter steget har tagits.
    ny_energi = energi - energiförlust

    # Om partikelns energi < 0, ansätt att partikelns energi = 0.
    if ny_energi <= 0:
        ny_energi = 0

    return ny_energi

"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FÖRFLYTTNING
Ta steg.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


@jit(nopython=True)
def förflyttning(x, y, z, dx, dy, dz, voxel_sidlängd=1):
    x_ny = x + dx / voxel_sidlängd
    y_ny = y + dy / voxel_sidlängd
    z_ny = z + dz / voxel_sidlängd

    x_round = round(x_ny)
    y_round = round(y_ny)
    z_round = round(z_ny)

    return x_ny, y_ny, z_ny, x_round, y_round, z_round


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FÖLJ PARTIKEL
Följer partikeln allteftersom den växelverkar.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def laddad_partikel_väg(energi_start, position_start, phi, theta, steglängd, radie_sfär, rho_medium,
                        stopping_power_data,
                        max_antal_steg):
    """
    Funktion som följer alfapartikeln allteftersom den växelverkar i ett medium.
    :param energi_start: Startenergi.
    :param position_start: Startposition.
    :param phi: Sfärisk vinkel.
    :param theta: Sfärisk vinkel.
    :param steglängd: Steglängd för partikeln i mediumet.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param rho_medium: Mediumets densitet.
    :param stopping_power_data: Tabellerad data för stopping power.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Energideponeringen innanför sfären.
    """

    """
    # Initiera tomma listor för att spara datan.
    x_list, y_list, z_list, dos_list = [], [], [], []

    position_vektor = position_start
    x, y, z = position_vektor[0], position_vektor[1], position_vektor[2]

    energi = energi_start

    # Dela upp steglängden i mindre steg.
    steg_storlek = steglängd / max_antal_steg

    # Håll reda på antalet steg.
    steg_tagna = 0

    # Under tiden som partikeln fortfarande inte tagit hela sitt steg,
    # samt att partikeln inte förlorat all sin energi.
    while steg_tagna < max_antal_steg and energi > 0:

        steg_tagna += 1

        # Ta ett litet steg.
        dx, dy, dz = transformera_koordinatsystem(steg_storlek, phi, theta, steg_storlek, 0, 0)
        x, y, z, _, _, _ = förflyttning(x, y, z, dx, dy, dz)
        position_vektor = np.array([x, y, z])

        # Håll reda på ifall partikeln befinner sig i sfären eller inte.
        # Ekvationen för en sfär: x^2 + y^2 + z^2 = r^2
        if np.sqrt(np.dot(position_vektor, position_vektor)) > radie_sfär:
            # print(f'Energideponering i position ', position_vektor)
            break

        else:
            # Beräkna energin efter steget.
            energi_ny = energi_efter_energiförlust(energi, steg_storlek, rho_medium, stopping_power_data)
            # print(f'ny energi: {energi * 10 ** (-6):.2f} MeV')

            # Spara mätpunkter.
            dos_list.append(energi - energi_ny)
            x_list.append(position_vektor[0])
            y_list.append(position_vektor[1])
            z_list.append(position_vektor[2])

            energi = energi_ny
            continue

    """
    # Initiera tomma listor för att spara datan.
    x_list, y_list, z_list, dos_list = [], [], [], []

    position_vektor = position_start
    energi = energi_start

    # Dela upp steglängden i mindre steg.
    steg_storlek = steglängd / max_antal_steg

    # Skapa en riktningsvektor.
    riktning = np.array(
        [np.sin(theta) * np.cos(phi)
            , np.sin(theta) * np.sin(phi)
            , np.cos(theta)])

    riktning /= np.linalg.norm(riktning)

    # Vektor för de små stegen som ska tas.
    steg_vektor = riktning * steg_storlek

    # Håll reda på antalet steg.
    steg_tagna = 0

    # Under tiden som partikeln fortfarande inte tagit hela sitt steg,
    # samt att energin är över tröskelenergin.
    while steg_tagna < max_antal_steg and energi > 1000:

        steg_tagna += 1

        # Ta ett litet steg.
        position_vektor += steg_vektor

        # Håll reda på ifall partikeln befinner sig i sfären eller inte.
        # Ekvationen för en sfär: x^2 + y^2 + z^2 = r^2
        if np.sqrt(np.dot(position_vektor, position_vektor)) > radie_sfär:
            # print(f'Energideponering i position ', position_vektor)
            break

        else:
            # Beräkna energin efter steget.
            energi_ny = energi_efter_energiförlust(energi, steg_storlek, rho_medium, stopping_power_data)
            # print(f'ny energi: {energi * 10 ** (-6):.2f} MeV')

            # Spara mätpunkter.
            dos_list.append(energi - energi_ny)
            x_list.append(position_vektor[0])
            y_list.append(position_vektor[1])
            z_list.append(position_vektor[2])

            energi = energi_ny
            continue


    # Beräkna den totala energideponeringen (i sfären) från partikeln.
    energideponering = energi_start - energi
    # print(f'energideponering: {energideponering} eV')
    return energideponering, x_list, y_list, z_list, dos_list


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIMULERING
Simuleringskod som sammanställer alla ovan nämnda avsnitt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""


def run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_alpha, radie_sfär,
                 max_antal_steg, fig_title):
    """
    Monte-Carlo simulering för alfapartiklarna.
    :param iterationer: Antal sönderfall som ska simuleras.
    :param rho_medium: Densiteten av mediumet (vatten).
    :param radie_partikel: Partikelns radie (alfa).
    :param stopping_power_data: Stopping power data.
    :param position_start_alpha: Fördelningsfunktion: uniform eller ytfördelning.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Summeringen av energideponeringen innanför sfären.
    """

    # Initiera en energisumma och tomma listor för att spara datan.
    energideponering_summa = 0
    x_list, y_list, z_list, dos_list = [], [], [], []

    #   -----------------------------------
    #   Vilken fördelningsfunktion som ska användas bestämmer hur
    #   sampling av riktning och position sker.
    #   -----------------------------------
    if position_start_alpha == position_start_skal:
        #   -----------------------------------
        #   Ytfördelning på en sfär.
        #   -----------------------------------
        iterationer = 0.5 * iterationer
        for i in range(int(iterationer)):
            # Sampla startenergin.
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6):.2f} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_skal()
            position_start = position_start_skal(radie_sfär, radie_partikel)

            # Sampla steglängd för partikeln.
            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering, x, y, z, dos = laddad_partikel_väg(energi, position_start, phi, theta, steglängd,
                                                                 radie_sfär,
                                                                 rho_medium, stopping_power_data, max_antal_steg)

            # Summera alla dosbidrag.
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6):.2f} MeV')

            # Spara mätpunkter för plottning.
            x_list += x
            y_list += y
            z_list += z
            dos_list += dos

    else:
        #   -----------------------------------
        #   Uniform fördelning i en sfär.
        #   -----------------------------------
        for i in range(iterationer):
            # Sampla startenergin.
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6):.2f} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_uniform()
            position_start = position_start_innanför(radie_sfär)

            # Sampla steglängd för partikeln.
            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering, x, y, z, dos = laddad_partikel_väg(energi, position_start, phi, theta, steglängd,
                                                                 radie_sfär,
                                                                 rho_medium, stopping_power_data, max_antal_steg)

            # Summera alla dosbidrag.
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6):.2f} MeV')

            # Spara mätpunkter för plottning.
            x_list += x
            y_list += y
            z_list += z
            dos_list += dos

    # print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')

    #   -----------------------------------
    #   Visualisera resultat i figur.
    #   -----------------------------------
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x_list, y_list, z_list, c=dos_list, cmap='plasma', label='Partikel position')
    # Fixa colorbar för att se energideponeringen i figuren

    # fig.colorbar(ax=ax, label='Energideponering',)

    ax.set_xlabel('x-axel (m)')
    ax.set_ylabel('y-axel (m)')
    ax.set_zlabel('z-axel (m)')

    # Testar att sätta en sfär för tumören
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = radie_sfär * np.cos(u) * np.sin(v)
    y = radie_sfär * np.sin(u) * np.sin(v)
    z = radie_sfär * np.cos(v)
    ax.plot_wireframe(x, y, z, color="k", alpha=0.3, label='Tumören')
    ax.legend(fontsize=font_size)
    plt.title(fig_title, fontsize=font_size_title)

    # Visa figur
    plt.tight_layout()
    plt.savefig(fig_title)
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
    #   -----------------------------------
    #   Kör simuleringen med ingångsvärden.
    #   -----------------------------------
    iterationer = 10 ** 2
    dummy_iterationer = 10 ** 1
    max_antal_steg = 10 ** 3

    fig_title_skal = 'Alfasönderfall vid ytfördelning'
    fig_title_innanför = 'Alfasönderfall vid uniform fördelning'

    stopping_power_data = np.loadtxt(stopping_power_alfa_file)

    rho_medium = rho_vatten
    radie_partikel = radie_alpha

    radie_sfär_skal = 300 * 10 ** (-6)
    radie_sfär_innanför = 1 * 10 ** (-3)

    print(
        '\n-----------------------------------\nDUMMY\n-----------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_skal,
                     radie_sfär_skal, max_antal_steg, fig_title_skal)

    start = time.time()

    print(
        '\n-----------------------------------\nRIKTIG\n-----------------------------------\n')
    energideponering_tot_skal = run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                             position_start_skal, radie_sfär_skal, max_antal_steg, fig_title_skal)
    energideponering_skal_Gy = energideponering_eV_till_Gy(energideponering_tot_skal, rho_medium, radie_sfär_skal)

    end_time(start)

    print(
        '\n-----------------------------------\nDUMMY\n-----------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_innanför,
                     radie_sfär_innanför, max_antal_steg, fig_title_innanför)

    start = time.time()
    print(
        '\n-----------------------------------\nRIKTIG\n-----------------------------------\n')
    energideponering_tot_innanför = run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                 position_start_innanför, radie_sfär_innanför, max_antal_steg,
                                                 fig_title_innanför)
    energideponering_innanför_Gy = energideponering_eV_till_Gy(energideponering_tot_innanför, rho_medium,
                                                               radie_sfär_innanför)
    end_time(start)

    print(
        '\n-----------------------------------\nRESULTAT\n-----------------------------------\n')

    #   -----------------------------------
    #   Beräkna resultat och jämför med valideringsdata.
    #   -----------------------------------
    print(
        f'\nSkal (300 mikrometer): Energideponering:\n{energideponering_skal_Gy * 10 ** 6 / iterationer} E-06 Gy / sönderfall')
    print(f'faktor {(energideponering_skal_Gy * 10 ** 6 / iterationer) / 1.66:.3f} av facit')

    print(
        f'\nInnanför (1 mm): Energideponering:\n{energideponering_innanför_Gy * 10 ** 8 / iterationer} E-08 Gy / sönderfall')
    print(f'faktor {(energideponering_innanför_Gy * 10 ** 8 / iterationer) / 9.18:.3f} av facit')
