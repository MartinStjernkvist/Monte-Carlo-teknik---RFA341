"""
----------------------------------------------------------------------
----------------------------------------------------------------------
IMPORTS
Pythonpaket och konstanter.
----------------------------------------------------------------------
----------------------------------------------------------------------
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

radie_alpha = 1.2 * 10 ** (-15) * 4 ** (1 / 3)  # Radie alfapartikel i meter (Physics Handbook)

At211_energi_MeV = [5.869, 5.2119, 5.1403, 4.9934, 4.895]  # Sönderfallsdata för At-211.
At211_energi = [(lambda x: x * 10 ** 6)(x) for x in At211_energi_MeV]  # Energi i eV
At211_intensitet = [41.78, 0.0039, 0.0011, 0.0004, 0.00004]  # Intensitet i % för energierna
At211_sannolikhet = np.cumsum(At211_intensitet) / np.sum(At211_intensitet)

rho_vatten = 998  # Vattens densitet i kg / m^3.

#   ----------------------------------------------------------------------
#   Filer med Data
#   ----------------------------------------------------------------------

stopping_power_alfa_file = 'Stoppingpower_data_alfa'


#   ----------------------------------------------------------------------
#   Funktioner
#   ----------------------------------------------------------------------

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
----------------------------------------------------------------------
----------------------------------------------------------------------
RIKTNING
Samplar riktningen för en partikel, utifrån fördelningsfunktion.
----------------------------------------------------------------------
----------------------------------------------------------------------
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
    phi = pi / 4 * (np.random.rand() + 2)
    return theta, phi


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
STARTPOSITION
Samplar startposition för partiklarna utifrån fördelningsfunktion.
----------------------------------------------------------------------
----------------------------------------------------------------------
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
----------------------------------------------------------------------
----------------------------------------------------------------------
STOPPING POWER OCH STEGLÄNGD
Beräkna stopping power och steglängden utifrån energi.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


def stopping_power_och_steglängd(energi, rho_medium, stopping_power_data):
    """
    Funktion som ger stopping power och
    :param energi: Partikelns energi i eV.
    :param rho_medium: Mediumets densitet i kg / m^3.
    :param stopping_power_data: Tabellerad data.
    :return: Stopping power (eV/m^2) och steglängden (m).
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
----------------------------------------------------------------------
----------------------------------------------------------------------
STARTPOSITION
Samplar startposition för partiklarna utifrån fördelningsfunktion.
----------------------------------------------------------------------
----------------------------------------------------------------------
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
----------------------------------------------------------------------
----------------------------------------------------------------------
STARTENERGI
Beräknar startenergin för partikeln, utifrån söndefallsdata.
----------------------------------------------------------------------
----------------------------------------------------------------------
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
----------------------------------------------------------------------
----------------------------------------------------------------------
ENERGIFÖRLUST
Beräknar energiförlusten efter varje steg som partikeln tar.
----------------------------------------------------------------------
----------------------------------------------------------------------
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
----------------------------------------------------------------------
----------------------------------------------------------------------
FÖLJ PARTIKEL
Följer partikeln allteftersom den växelverkar.
----------------------------------------------------------------------
----------------------------------------------------------------------
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

    position_vektor = position_start
    energi = energi_start

    # Dela upp steglängden i mindre steg.
    steg_storlek = steglängd / max_antal_steg

    # Skapa en riktningsvektor.
    riktning = np.array(
        [np.sin(theta) * np.cos(phi)
            , np.sin(theta) * np.cos(phi)
            , np.cos(theta)])

    riktning /= np.linalg.norm(riktning)

    # Vektor för de små stegen som ska tas.
    steg_vektor = riktning * steg_storlek

    # Håll reda på antalet steg.
    steg_tagna = 0

    # Under tiden som partikeln fortfarande inte tagit hela sitt steg,
    # samt att partikeln inte förlorat all sin energi.
    while steg_tagna < max_antal_steg and energi > 0:

        steg_tagna += 1

        # Ta ett litet steg.
        position_vektor += steg_vektor

        # Beräkna energin efter steget.
        energi = energi_efter_energiförlust(energi, steg_storlek, rho_medium, stopping_power_data)
        # print(f'ny energi: {energi * 10 ** (-6):.2f} MeV')

        # Håll reda på ifall partikeln befinner sig i sfären eller inte.
        # Ekvationen för en sfär: x^2 + y^2 + z^2 = r^2
        if np.sqrt(np.dot(position_vektor, position_vektor)) < radie_sfär:
            # print(f'Energideponering i position ', position_vektor)
            continue
        else:
            break
        # if np.dot(position_vektor, position_vektor) > radie_sfär:
        #     break

    # Beräkna den totala energideponeringen (i sfären) från partikeln.
    energideponering = energi_start - energi
    # print(f'energideponering: {energideponering} eV')
    return energideponering


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
MONTE-CARLO SIMULERING
Samplar startposition för partiklarna utifrån fördelningsfunktion.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


def run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_alpha, radie_sfär,
                 max_antal_steg):
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

    energideponering_summa = 0
    #   ----------------------------------------------------------------------
    #   Vilken fördelningsfunktion som ska användas bestämmer hur
    #   sampling av riktning och position sker.
    #   ----------------------------------------------------------------------
    if position_start_alpha == position_start_skal:
        #   ----------------------------------------------------------------------
        #   Ytfördelning på en sfär.
        #   ----------------------------------------------------------------------
        iterationer = 0.5 * iterationer
        for i in range(int(iterationer)):
            # Sampla startenergin.
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6)} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_skal()
            position_start = position_start_skal(radie_sfär, radie_partikel)

            # Sampla steglängd för partikeln.
            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering = laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, stopping_power_data, max_antal_steg)

            # Summera alla dosbidrag.
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    else:
        #   ----------------------------------------------------------------------
        #   Uniform fördelning i en sfär.
        #   ----------------------------------------------------------------------
        for i in range(iterationer):
            # Sampla startenergin.
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6)} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_uniform()
            position_start = position_start_innanför(radie_sfär)

            # Sampla steglängd för partikeln.
            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering = laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, stopping_power_data, max_antal_steg)

            # Summera alla dosbidrag.
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    # print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
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
    #   Kör simuleringen med ingångsvärden.
    #   ----------------------------------------------------------------------
    iterationer = 10 ** 2
    dummy_iterationer = 10 ** 1
    max_antal_steg = 10 ** 3

    stopping_power_data = np.loadtxt(stopping_power_alfa_file)

    rho_medium = rho_vatten
    radie_partikel = radie_alpha

    radie_sfär_skal = 300 * 10 ** (-6)
    radie_sfär_innanför = 1 * 10 ** (-3)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_skal,
                     radie_sfär_skal, max_antal_steg)

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                             position_start_skal, radie_sfär_skal, max_antal_steg)
    energideponering_skal_Gy = energideponering_eV_till_Gy(energideponering_tot_skal, rho_medium, radie_sfär_skal)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_innanför,
                     radie_sfär_innanför, max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                 position_start_innanför, radie_sfär_innanför, max_antal_steg)
    energideponering_innanför_Gy = energideponering_eV_till_Gy(energideponering_tot_innanför, rho_medium,
                                                               radie_sfär_innanför)
    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nRESULTAT\n----------------------------------------------------------------------\n')

    #   ----------------------------------------------------------------------
    #   Beräkna resultat och jämför med valideringsdata.
    #   ----------------------------------------------------------------------
    print(
        f'\nSkal (300 mikrometer): Energideponering:\n{energideponering_skal_Gy * 10 ** 6 / iterationer} E-06 Gy / sönderfall')
    print(f'faktor {(energideponering_skal_Gy * 10 ** 6 / iterationer) / 1.66} av facit')

    print(
        f'\nInnanför (1 mm): Energideponering:\n{energideponering_innanför_Gy * 10 ** 8 / iterationer} E-08 Gy / sönderfall')
    print(f'faktor {(energideponering_innanför_Gy * 10 ** 8 / iterationer) / 9.18} av facit')
