"""
----------------------------------------------------------------------
----------------------------------------------------------------------
IMPORTS
Pythonpaket, konstanter, och några generella funktioner.
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

radie_alpha = 1.9 * 10 ** (-15)  # m

massa_H = 1  # enhet u
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

tvärsnitt_file = '../given_data/Tvärsnittstabeller_Fotoner.xlsx'
attenueringsdata_file = '../given_data/Attenueringsdata.xlsx'
anatomidefinitioner_file = '../given_data/Anatomidefinitioner.xlsx'

mat_file = 'phantom_data.mat'


# tvärsnitt_file = r'C:/Users/Admin/Documents/GitHub/Monte Carlo Linnea/Tvärsnittstabeller_Fotoner.xlsx'
# attenueringsdata_file = r"C:/Users/Admin/Documents/GitHub/Monte Carlo Linnea/Attenueringsdata.xlsx"
# anatomidefinitioner_file = r"C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea/Anatomidefinitioner.xlsx"
# mat_file = r"C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\phantom_data.mat"


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
----------------------------------------------------------------------
----------------------------------------------------------------------
VISUALISERA BIN FIL
Visualisera fantommatrisen.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


class visualisera_matris:
    """
    Kod som visualiserar fantomen i tre olika tvärsnitt, med sliders.
    """

    def __init__(self, array_3d, visa_något=False):
        self.array_3d = array_3d
        self.x, self.y, self.z = array_3d.shape

        # Initialisera mitt-tvärsnittet för respektive dimension.
        self.slice_x_index = self.x // 2
        self.slice_y_index = self.y // 2
        self.slice_z_index = self.z // 2

        # Skapa figur.
        self.fig = plt.figure(figsize=(15, 7.5))
        self.gs = GridSpec(1, 3, width_ratios=[1, 1, 1], height_ratios=[1])

        # Skapa delfigurer.
        self.ax1 = plt.subplot(self.gs[0])
        self.ax2 = plt.subplot(self.gs[1])
        self.ax3 = plt.subplot(self.gs[2])

        # visa_något bestämmer ifall det som plottas ska ha olika färger beroende på voxelvärde, eller
        # ifall allt som plottas ska vara vitt (så att det syns)
        if visa_något == False:
            self.vmin, self.vmax = 0, np.max(array_3d)
        else:
            self.vmin, self.vmax = 0, 0.0001 * np.max(array_3d)

        self.img1 = self.ax1.imshow(array_3d[self.slice_x_index, :, :], cmap='gray', vmin=self.vmin,
                                    vmax=self.vmax)  # x-tvärsnitt
        self.ax1.set_title(f'X Slice: {self.slice_x_index}')
        self.img2 = self.ax2.imshow(array_3d[:, self.slice_y_index, :], cmap='gray', vmin=self.vmin,
                                    vmax=self.vmax)  # y-tvärsnitt
        self.ax2.set_title(f'Y Slice: {self.slice_y_index}')
        self.img3 = self.ax3.imshow(array_3d[:, :, self.slice_z_index], cmap='gray', vmin=self.vmin,
                                    vmax=self.vmax)  # z-tvärsnitt
        self.ax3.set_title(f'Z Slice: {self.slice_z_index}')

        # Skapa sliders.
        self.create_sliders()

        # Attach update function to sliders
        self.slider_x.on_changed(self.update)
        self.slider_y.on_changed(self.update)
        self.slider_z.on_changed(self.update)

    def create_sliders(self):
        ax_slider_x = plt.axes([0.1, 0.01, 0.75, 0.03], facecolor='lightgoldenrodyellow')
        self.slider_x = Slider(ax_slider_x, 'X Slice', 0, self.x - 1, valinit=self.slice_x_index, valstep=1)

        ax_slider_y = plt.axes([0.1, 0.06, 0.75, 0.03], facecolor='lightgoldenrodyellow')
        self.slider_y = Slider(ax_slider_y, 'Y Slice', 0, self.y - 1, valinit=self.slice_y_index, valstep=1)

        ax_slider_z = plt.axes([0.1, 0.11, 0.75, 0.03], facecolor='lightgoldenrodyellow')
        self.slider_z = Slider(ax_slider_z, 'Z Slice', 0, self.z - 1, valinit=self.slice_z_index, valstep=1)

    def update(self, val):
        # Updatera sliders.
        self.slice_x_index = int(self.slider_x.val)
        self.slice_y_index = int(self.slider_y.val)
        self.slice_z_index = int(self.slider_z.val)

        # Updatera figurerna.
        self.img1.set_data(self.array_3d[self.slice_x_index, :, :])
        self.ax1.set_title(f'X Slice: {self.slice_x_index}')

        self.img2.set_data(self.array_3d[:, self.slice_y_index, :])
        self.ax2.set_title(f'Y Slice: {self.slice_y_index}')

        self.img3.set_data(self.array_3d[:, :, self.slice_z_index])
        self.ax3.set_title(f'Z Slice: {self.slice_z_index}')

        self.fig.canvas.draw_idle()

    def show(self):
        plt.show()


#   ----------------------------------------------------------------------
#   Läs in data
#   ----------------------------------------------------------------------

"""
I matlab:

% Save the entire array_3d to a .mat file
save('phantom_data.mat', 'array_3d');
"""

data = scipy.io.loadmat(mat_file)
fantom_matris = data['array_3d']

"""
----------------------------------------------------------------------
----------------------------------------------------------------------
MATRISER
Kapa fantommatrisen, skapa en njurmatris och benmärgmatris.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""

#   ----------------------------------------------------------------------
#   Skapa kapade matriser, utifrån fantommatrisen.
#   ----------------------------------------------------------------------

slicad_fantom_matris = fantom_matris[:, 34:212, 450:1182]  # inkluderar fantom, avskuret vid benen

x, y, z = slicad_fantom_matris.shape

#   ----------------------------------------------------------------------
#   Njure.
#   ----------------------------------------------------------------------

# Skapa tom matris för att sedan fylla i med värden.
slicad_njure_matris = np.zeros((x, y, z))

# Genom att ansätta värden för njur-vävnad från filen Anatomidefinitioner.xlsx kan en matris skapas
# som filtrerar bort alla voxelvärden från fantommatrisen som inte motsvarar njurvävnad.
target_values = [17, 18, 19, 20, 21, 22, 23]
mask = np.isin(slicad_fantom_matris, target_values)
slicad_njure_matris[mask] = slicad_fantom_matris[mask]

#   ----------------------------------------------------------------------
#   Benmärg.
#   ----------------------------------------------------------------------
"""
Samma grej som för njurarna, skillnaden är att vi byter target_values till värden med benmärg
"""

slicad_benmärg_matris = np.zeros((x, y, z))

# Samma grej som för njurarna, skillnaden är att vi byter target_values till värden med benmärg
target_values = [29]
mask = np.isin(slicad_fantom_matris, target_values)
slicad_benmärg_matris[mask] = slicad_fantom_matris[mask]

"""
----------------------------------------------------------------------
----------------------------------------------------------------------
SAMPLA ENERGI
Sampla startenergin för fotonerna.
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
SAMPLA STARTPOSITION
Sampla startposition för fotonerna (i njuren).
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


@jit(nopython=True)
def position_start(slicad_njure_matris):
    """
    Funktion som samplar en startposition utifrån en matris med endast njurvävnad.
    :param slicad_njure_matris: Matris med njurvävnad.
    :return: Slumpad startposition
    """
    x_size, y_size, z_size = slicad_njure_matris.shape

    # Tar första möjliga voxel vars värde inte är 0.
    while True:
        x = np.random.randint(0, x_size)
        y = np.random.randint(0, y_size)
        z = np.random.randint(0, z_size)

        if slicad_njure_matris[x, y, z] != 0:
            break

    return x, y, z  # , x_round, y_round, z_round


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
SAMPLA RIKTNING OCH STEG
Sampla riktning och ta steg i en viss riktning.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


@jit(nopython=True)
def riktning_uniform():
    """
    Funktion som samplar uniformt spridningsvinklarna för en foton.
    Används både till startpositionen, samt vid fotoabsorption.
    :return: Spridningsvinklar theta och phi (sfäriska koordinater).
    """

    # cos(theta) ska vara mellan -1 och 1
    theta = np.arccos(-1 + 2 * np.random.rand())

    phi = 2 * pi * np.random.rand()
    return theta, phi


@jit(nopython=True)
def steg(theta, phi, steglängd):
    """
    Funktion som tar ett steg i en specificerad riktning.
    :param theta: Spridningsvinkel.
    :param phi: Spridningsvinkel.
    :param x: Startposition x.
    :param y: Startposition y.
    :param z: Startposition z.
    :return: Ny position.
    """

    # Steg, i termer av voxelsidlängd.
    dx = steglängd * np.sin(theta) * np.cos(phi) / voxel_sidlängd
    dy = steglängd * np.sin(theta) * np.sin(phi) / voxel_sidlängd
    dz = steglängd * np.cos(theta) / voxel_sidlängd

    return dx, dy, dz


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
SAMPLA STEGLÄNGD
Sampla steglängd.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


@jit(nopython=True)
def invers_funktion(x, mu):
    """
    Invers funktion, som används för att sampla fotonens steglängd.
    :param x: Ett slumpat tal mellan 0 och 1.
    :param mu: Attenueringskoefficient
    """
    return -np.log(x) / mu


@jit(nopython=True)
def medelvägslängd(mu):
    """
    Funktion som samplar steglängden utifrån den inverstransformerade funktionen ovan.
    """
    medelvägslängd = invers_funktion(np.random.rand(), mu) / voxel_sidlängd  # LÄGG TILL VOXELLÄNGD
    return medelvägslängd


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
ATTENUERINGSDATA
Bestäm attenueringskoefficienterna.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


class attenueringsdata:

    def __init__(self, voxelvärde, energi, df_attenueringsdata, df_anatomidefinitioner):

        # Läs in data från argumenten till klassen. Excelfilerna bearbetas med paketet pandas.
        self.df_attenueringsdata = df_attenueringsdata
        self.anatomidefinitioner = df_anatomidefinitioner

        self.voxelvärde = voxelvärde
        self.energi = energi

        self.energi_list = self.df_attenueringsdata['E'].to_list()

    def voxelvärde_till_material(self):
        """
        Funktion som omvandlar voxelvärdet i fantommatrisen till ett material.
        Voxelvärdet motsvarar viss sorts vävnad, enligt filen Anatomidefinitioner.xlsx.
        Därför behövs en koppling mellan vävnaderna i anatomidefinitioner och materialen i filen Attenueringsdata.xlsx.
        :return: 1) materialnamnet som motsvarar ett visst voxelvärde och 2) en namnlista med material som ingår i fantomen.
        """
        # print(vävnad)

        # luft = 0
        # bihåla? antag vatten
        # matstrupe? antag vatten
        # svalg? antag muskler
        # kortikalt ben? antag röd benmärg
        # benmärg? antingen röd eller yellow, antag röd
        # prostata? antar bladder
        # skull? antag benmärg, antag röd
        # bronker? antag vatten

        # Skapa listor för materialen, där ett material motsvarar voxelvärden (vävnad enligt Anatomidefinitioner).
        water = [3, 8, 12, 37]
        muscle = [6, 13]
        lung = [11]
        dry_spine = [27, 28]
        dry_rib = [25]
        # adipose = [] # Hittade ingenting som passade in här
        blood = [2]
        heart = [1]
        kidney = [17, 18, 19, 20, 21, 22, 23]
        liver = [9]
        lymph = [36, 41]
        pancreas = [16]
        intestine = [14, 15, 32, 33]
        # skull = [] # Hittade ingenting som passade in här
        cartilage = [34]
        brain = [7]
        spleen = [24]
        air = [0, 35, 38]
        breast_mammary = [5]
        skin = [4]
        eye_lens = [42, 43]
        # ovary = [] # Hittade ingenting som passade in här
        red_marrow = [26, 29]
        yellow_marrow = [44]
        # testis = [] # Hittade ingenting som passade in här
        thyroid = [39, 40]
        bladder = [10, 30, 31]

        # Sätt ihop individuella materiallistor till en stor lista.
        big_list = [water, muscle, lung, dry_spine, dry_rib, blood, heart, kidney, liver, lymph, pancreas, intestine,
                    cartilage, brain, spleen, air, breast_mammary, skin, eye_lens, red_marrow, yellow_marrow, thyroid,
                    bladder]

        # Lista med namn, för att kunna indexera i attenueringsdatan.
        big_list_names = ['water', 'muscle', 'lung', 'dry_spine', 'dry_rib', 'blood', 'heart', 'kidney', 'liver',
                          'lymph', 'pancreas', 'intestine',
                          'cartilage', 'brain', 'spleen', 'air', 'breast_mammary', 'skin', 'eye_lens', 'red_marrow',
                          'yellow_marrow', 'thyroid',
                          'bladder']

        # Filtrerar en materiallista i taget i den större, sammansatta listan.
        # Om voxelvärdet matchar ett av värdena i en materiallista erhålls materialnamnet.
        for i in range(len(big_list)):
            sublist = big_list[i]
            if any(self.voxelvärde == värde for värde in sublist):
                material = big_list_names[i]
                break

        return material, big_list_names

    def mu(self):
        """
        Funktion som tar fram attenueringskoefficienten utifrån energin och voxelvärdet i en position.
        :return: attenueringskoefficient
        """
        # Energilistan (första kolumnen i Attenueringsdata.xlsx) läses in.
        # De två närmaste indexen, för energierna i listan som är närmast den ingående energin för funktionen tas fram.
        energi_list = np.array(self.energi_list)
        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]

        # Kalla på ovanstående funktion, för att omvandla voxelvärde till material.
        material, _ = self.voxelvärde_till_material()
        mu_list = np.array(self.df_attenueringsdata[material].to_list())

        # Närmaste energierna och attenueringskoefficienterna.
        energi_close = energi_list[closest_indices]
        mu_close = mu_list[closest_indices]

        # Om närliggande attenueringskoefficienter är lika stora -> attenueringskoefficient = första värdet.
        if mu_close[1] - mu_close[0] < 10 ** (-15):
            mu_target = mu_close[0]

        # Om närliggande attenueringskoefficienter inte är lika stora -> linjärinterpolera fram attenueringskoefficienten.
        else:
            mu_target = mu_close[0] + (self.energi - energi_close[0]) * (mu_close[1] - mu_close[0]) / (
                    energi_close[1] - energi_close[0])

        # print(f'mu target {mu_target}')

        return mu_target

    def mu_max(self):
        """
        Funktion som tar fram den maximala attenueringskoefficienten i fantomen för en viss energi.
        :return: maximal attenueringskoefficient
        """
        # Energilistan (första kolumnen i Attenueringsdata.xlsx) läses in.
        # De två närmaste indexen, för energierna i listan som är närmast den ingående energin för funktionen tas fram.
        energi_list = np.array(self.energi_list)
        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]

        # Kalla på ovanstående funktion, för att erhålla namnlistan med material som ingår i fantomen.
        _, big_list_names = self.voxelvärde_till_material()

        # Skapa en tom vektor, som sedan ska fyllas.
        mu_array = np.zeros(len(big_list_names), dtype=object)

        # För varje material i namnlistan tas värden på attenueringskoefficienten för samma "närmaste index" som för energin.
        # En vektor med listor skapas.
        for i in range(len(big_list_names)):
            mu_array[i] = np.array(self.df_attenueringsdata[big_list_names[i]][closest_indices].to_list())

        # Sortera bland vektorn med listorna, för att ta fram materialet (indexet) som har störst attenueringskoefficient.
        mu_max_close_values = [np.max(arr) for arr in mu_array]
        mu_max_close = max(mu_max_close_values)
        mu_max_close_index = mu_max_close_values.index(mu_max_close)

        # Ta ut listan med de två attenueringskoefficienterna, för materialet med högst attenueringskoefficient.
        mu_max_close = mu_array[mu_max_close_index]
        # Även energierna som ligger närmast energin som matas in i funktionen.
        energi_close = energi_list[closest_indices]

        # Om närliggande attenueringskoefficienter är lika stora -> attenueringskoefficient = första värdet.
        if mu_max_close[1] - mu_max_close[0] < 10 ** (-15):
            mu_max = mu_max_close[0]

        # Om närliggande attenueringskoefficienter inte är lika stora -> linjärinterpolera fram attenueringskoefficienten.
        else:
            mu_max = mu_max_close[0] + (self.energi - energi_close[0]) * (mu_max_close[1] - mu_max_close[0]) / (
                    energi_close[1] - energi_close[0])

        return mu_max


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
VÄXELVERKAN
Bestäm vilken växelverkanprocess som sker för en viss energi.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


class växelverkan:

    def __init__(self, energi, df_tvärsnitt):
        self.df = df_tvärsnitt

        self.energi = energi

        self.energi_list = self.df['Energy (eV)'].to_list()
        self.foto_list = self.df['Photoelectric  (cm^2)'].to_list()
        self.compton_list = self.df['Compton (cm^2)'].to_list()
        self.rayleigh_list = self.df['Rayleigh (cm^2)'].to_list()

    def find_foto_tvärsnitt(self):
        energi_list = np.array(self.energi_list)
        foto_list = np.array(self.foto_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        foto_close = foto_list[closest_indices]
        energi_close = energi_list[closest_indices]

        if energi_close[1] - energi_close[0] < 10 ** (-15):
            foto_target = foto_close[0]

        else:
            # linjär interpolering funktion
            foto_target = foto_close[0] + (self.energi - energi_close[0]) * (foto_close[1] - foto_close[0]) / (
                    energi_close[1] - energi_close[0])
        return foto_target

    def find_compton_tvärsnitt(self):
        energi_list = np.array(self.energi_list)
        compton_list = np.array(self.compton_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        compton_close = compton_list[closest_indices]
        energi_close = energi_list[closest_indices]

        if energi_close[1] - energi_close[0] < 10 ** (-15):
            compton_target = compton_close[0]

        else:
            compton_target = compton_close[0] + (self.energi - energi_close[0]) * (
                    compton_close[1] - compton_close[0]) / (
                                     energi_close[1] - energi_close[0])

        return compton_target

    def find_rayleigh_tvärsnitt(self):
        energi_list = np.array(self.energi_list)
        rayleigh_list = np.array(self.rayleigh_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        rayleigh_close = rayleigh_list[closest_indices]
        energi_close = energi_list[closest_indices]

        if energi_close[1] - energi_close[0] < 10 ** (-15):
            rayleigh_target = rayleigh_close[0]

        else:
            rayleigh_target = rayleigh_close[0] + (self.energi - energi_close[0]) * (
                    rayleigh_close[1] - rayleigh_close[0]) / (
                                      energi_close[1] - energi_close[0])
        return rayleigh_target

    def bestäm_växelverkan(self):
        foto_target = self.find_foto_tvärsnitt()
        compton_target = self.find_compton_tvärsnitt()
        rayleigh_target = self.find_rayleigh_tvärsnitt()

        tvärsnitt_lista = [foto_target, compton_target, rayleigh_target]  # cm^2

        tvärsnitt_lista_norm = np.cumsum(tvärsnitt_lista) / np.sum(tvärsnitt_lista)

        slump_tal = np.random.rand()
        if slump_tal <= tvärsnitt_lista_norm[0]:
            text = 'foto'
        elif slump_tal <= tvärsnitt_lista_norm[1]:
            text = 'compton'
        else:
            text = 'rayleigh'

        return text


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
FOTOABSORPTION
Sampla utfallet av fotoabsorption.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


@jit(nopython=True)
def foto_vxv(foton_energi):
    """
    Funktion som samplar utfallet av fotoväxelverkan.
    Antingen sker fluorescens eller inte, Augerelektroner antas deponera energi lokalt.
    :return: Energideponeringen som konsekvens av växelverkan, samt ifall fotonen slutar existera.
    """

    # Fluorescens.
    if np.random.rand() < fluorescence_yield:
        # X-ray
        energi_deponering = foton_energi - K_alpha
        attenuerad = 0
    else:
        # Augerelektron, lokal energideponering -> finns ingen foton att följa.
        energi_deponering = foton_energi
        attenuerad = 1

    return energi_deponering, attenuerad


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
COMPTONSPRIDNING
Sampla spridningsvinkel och energideponering.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


@jit(nopython=True)
def compton_vinkel_och_energiförlust(foton_energi):
    """
    Khans rejektionsalgoritm, artikeln "A Monte Carlo program for the simulation of scintillation camera characteristics"
    Algoritmen samplar vinkeln och energideponeringen för en foton som Comptonväxelverkar.
    :return: Spridningsvinkeln theta, den spridda fotonens energi och den spridda fotonens energiförlust (=energideponering).
    """

    while True:
        # Slumpa tre tal.
        R_1 = np.random.rand()
        R_2 = np.random.rand()
        R_3 = np.random.rand()

        # Termen alpha är fotonens energi i termer av vilomassan (energi) för en elektron.
        alpha = foton_energi / E_e

        if R_1 <= (2 * alpha + 1) / (2 * alpha + 9):

            eta = 2 * alpha * R_2
            if R_3 <= 4 * (eta ** -1 - eta ** -2):
                theta = np.arccos(1 - 2 * R_2)
                spridd_foton_energi = foton_energi / eta
                energi_förlust = foton_energi - spridd_foton_energi

                break

        elif R_1 > (2 * alpha + 1) / (2 * alpha + 9):

            eta = (2 * alpha + 1) / (2 * R_2 * alpha + 1)
            theta = np.arccos(1 - (eta - 1) / alpha)

            if R_3 <= 0.5 * (np.cos(theta) ** 2 + eta ** -1):
                spridd_foton_energi = foton_energi / eta
                theta = np.arccos(1 - (eta - 1) / alpha)
                energi_förlust = foton_energi - spridd_foton_energi

                break

        else:
            continue

    return theta, spridd_foton_energi, energi_förlust


"""
----------------------------------------------------------------------
----------------------------------------------------------------------


----------------------------------------------------------------------
----------------------------------------------------------------------
"""