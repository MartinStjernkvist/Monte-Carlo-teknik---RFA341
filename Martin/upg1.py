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
FÖRFLYTTNING
Ta ett steg, samt erhåll voxelkoordinater.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


@jit(nopython=True)
def förflyttning(x, y, z, dx, dy, dz):
    x = x + dx
    y = y + dy
    z = z + dz

    x_round = round(x)
    y_round = round(y)
    z_round = round(z)

    return x, y, z, x_round, y_round, z_round


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
TRANSFORMATION - HOMOGEN MATRIS
Translation och rotation i tre dimensioner.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


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

    # för att x-axeln ska sammanfalla med riktningsvektorn måste rotationsvinkeln vara pi/2 - theta
    angle = -(pi / 2 - theta_A)
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
    vektor_voxel = (1 / voxel_sidlängd) * vektor
    dx, dy, dz = vektor_voxel[0], vektor_voxel[1], vektor_voxel[2]

    return dx, dy, dz


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
BESTÄM OM INNANFÖR MATRIS
Bestäm om fotonen ska fortsätta följas eller inte.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


@jit(nopython=True)
def bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size, utanför_fantom, slicad_fantom_matris,
                         foton_energi):
    """
    Funktion som bestämmer om en foton ska fortsätta följas eller inte.
    Termen "attenuerad" i koden kan även vara att fotonen rymt från fantommatrisen.
    :param x_round: voxelposition x
    :param y_round: voxelposition y
    :param z_round: voxelposition z
    :param x_size: storlek matris x
    :param y_size: storlek matris y
    :param z_size: storlek matris z
    :param utanför_fantom: värde som ökar med 1, för att hålla reda på vad som händer när koden körs
    :param slicad_fantom_matris: fantommatrisen
    :return: om attenuerad = 1 kommer fotonen att sluta följas
    """

    # Om positionen för fotonen är utanför fantommatrisen, indikera att fotonen ska sluta följas.
    if (
            x_round < 0
            or x_round >= x_size
            or y_round < 0
            or y_round >= y_size
            or z_round < 0
            or z_round >= z_size
    ):
        utanför_fantom += 1
        attenuerad = 1

    # Om voxelvärdet för positionen är innanför matrisen, men befinner sig i luft (voxelvärde = 0)
    # eller
    # fotonenergin är mindre än gränsenergin för fotoabsorption
    # -> indikera att fotonen ska sluta följas.
    elif slicad_fantom_matris[x_round, y_round, z_round] == 0 or foton_energi < foton_energi_threshhold:
        utanför_fantom += 1
        attenuerad = 1

    else:
        attenuerad = 0

    return attenuerad, utanför_fantom


"""
----------------------------------------------------------------------
----------------------------------------------------------------------
BESTÄM OM VÄXELVERKAN
Bestäm om växelverkan sker vid en position, efter ett steg.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


def bestäm_om_vxv(voxelvärde, energi, mu_max, df_attenueringsdata, df_anatomidefinitioner):
    instans = attenueringsdata(voxelvärde, energi, df_attenueringsdata, df_anatomidefinitioner)
    mu = instans.mu()

    if np.random.rand() <= mu / mu_max:
        vxv_sker = True
    else:
        vxv_sker = False

    return vxv_sker


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
SIMULERING
Simuleringskod som sammanställer alla ovan nämnda avsnitt.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


def run_MC_multiprocess(args):
    (start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
     slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi, radionuklid_intensitet,
     radionuklid_sannolikhet) = args

    #   ----------------------------------------------------------------------
    #   Läs in data
    #   ----------------------------------------------------------------------
    df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)
    df_anatomidefinitioner = pd.read_excel(anatomidefinitioner_file, index_col=None)
    df_tvärsnitt = pd.read_excel(tvärsnitt_file, index_col=None)

    # Matrisdimensionerna, viktigt för att bedöma ifall fotonen är innanför matrisen eller inte.
    x_size, y_size, z_size = slicad_fantom_matris.shape

    # Skapa tom matris, för att sedan fylla på med energideponering i voxklarna.
    benmärg_matris_deponerad_energi = np.zeros((x_size, y_size, z_size))

    #   ----------------------------------------------------------------------
    #   Håll reda på vad som händer när koden körs
    #   ----------------------------------------------------------------------
    utanför_fantom = 0
    vxv_foto = 0
    träff_benmärg = 0

    #   ----------------------------------------------------------------------
    #   Loopar igenom alla iterationer
    #   ----------------------------------------------------------------------
    for i in range(start, end):

        # Initiera loopen med att attenuerad = 0.
        attenuerad = 0
        # print(i)

        # Start: sampla position, riktning och energi
        foton_energi = energi_start(radionuklid_energi, radionuklid_sannolikhet)
        x_start, y_start, z_start = position_start(slicad_njure_matris)
        theta, phi = riktning_uniform()

        voxelvärde = slicad_fantom_matris[x_start, y_start, z_start]
        instans = attenueringsdata(voxelvärde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
        mu_max = instans.mu_max()

        vxv_sker = False
        while vxv_sker == False:
            # Sampla medelvägslängden från inverstransformerad attenueringsfunktion.
            steglängd = medelvägslängd(mu_max)
            # print(f'steglängd: {steglängd}')

            # Gå steget till ny position från startpositionen i startriktningen.
            # Eftersom voxlarna har diskreta positioner måste avrundning till närmaste heltal göras.
            dx, dy, dz = steg(theta, phi, steglängd)
            x, y, z, x_round, y_round, z_round = förflyttning(x_start, y_start, z_start, dx, dy, dz)

            #   ----------------------------------------------------------------------
            #   Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen
            #   Ekvationen förekommer nedan efter varje nytt steg tas, dock utan if-statement.
            #   ----------------------------------------------------------------------

            attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size,
                                                              utanför_fantom, slicad_fantom_matris, foton_energi)

            if attenuerad == 0:
                voxelvärde = slicad_fantom_matris[x_round, y_round, z_round]

                vxv_sker = bestäm_om_vxv(voxelvärde, foton_energi, mu_max, df_attenueringsdata,
                                         df_anatomidefinitioner)
            else:
                # utanför matris
                vxv_sker = True

        # Kanske kolla om den fortsätter eller vxv?

        if attenuerad == 1:
            i += 1

        else:

            #   ----------------------------------------------------------------------
            #   Loopa under tiden som fotonen inte attenuerats och fortfarande är i matrisen.
            #   ----------------------------------------------------------------------
            while attenuerad == 0:
                # print('steglängd: ', steglängd)

                # Identifiera vilken voxel fotonen befinner sig i.
                voxelvärde = slicad_fantom_matris[x_round, y_round, z_round]
                instans = attenueringsdata(voxelvärde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
                mu_max = instans.mu_max()

                # Bestäm vilken typ av växelverkan som sker vid nya positionen.
                instans = växelverkan(foton_energi, df_tvärsnitt)
                vxv = instans.bestäm_växelverkan()
                # print(f'energi: {foton_energi * 10 ** (-3)} keV, vxv: {vxv}')

                #   ----------------------------------------------------------------------
                #   Fotoabsorption.
                #   ----------------------------------------------------------------------
                if vxv == 'foto':
                    vxv_foto += 1

                    # Bestäm om fluorescens sker eller inte.
                    energi_deponering, attenuerad = foto_vxv(foton_energi)
                    foton_energi = foton_energi - energi_deponering

                    # Registrera energideponering ifall fotonen växelverkar i en voxel med benmärg.
                    if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                        träff_benmärg += 1

                        benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energi_deponering

                        print(
                            f'foto: {energi_deponering:.0f} eV i voxel [{x_round, y_round, z_round}]')

                    # Om fluorescens sker -> följ ny foton.
                    if attenuerad == 0:
                        # Sampla spridningsvinklar (uniformt samplade).
                        theta_foto, phi_foto = riktning_uniform()

                        vxv_sker = False
                        while vxv_sker == False:
                            # Sampla medelvägslängden från inverstransformerad attenueringsfunktion.
                            steglängd_foto = medelvägslängd(mu_max)

                            dx, dy, dz = steg(theta, phi, steglängd_foto)
                            x, y, z, x_round, y_round, z_round = förflyttning(x, y, z, dx, dy, dz)

                            # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                            attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size,
                                                                              z_size, utanför_fantom,
                                                                              slicad_fantom_matris,
                                                                              foton_energi)

                            if attenuerad == 0:
                                voxelvärde = slicad_fantom_matris[x_round, y_round, z_round]

                                vxv_sker = bestäm_om_vxv(voxelvärde, foton_energi, mu_max, df_attenueringsdata,
                                                         df_anatomidefinitioner)
                            else:
                                vxv_sker = True

                        # Ingångsvärden till koordinat-transformeringen som behöver genomföras ifall växelverkan = compton eller rayleigh.
                        theta, phi = theta_foto, phi_foto
                        steglängd = steglängd_foto

                #   ----------------------------------------------------------------------
                #   Comptonspridning.
                #   ----------------------------------------------------------------------
                elif vxv == 'compton':

                    # Sampla spridningsvinklar, samt energideponeringen vid Comptonväxelverkan.
                    theta_compton, foton_energi, energideponering_compton = compton_vinkel_och_energiförlust(
                        foton_energi)
                    phi_compton = 2 * pi * np.random.rand()

                    # Registrera energideponering ifall fotonen växelverkar i en voxel med benmärg.
                    if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                        träff_benmärg += 1

                        benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energideponering_compton

                        print(
                            f'compton: {energideponering_compton:.0f} eV i voxel [{x_round, y_round, z_round}]')

                    vxv_sker = False
                    while vxv_sker == False:
                        # Sampla medelvägslängden från inverstransformerad attenueringsfunktion.
                        steglängd_compton = medelvägslängd(mu_max)

                        # Koordinattransformation, eftersom spridningsvinklarna för Comptonspridning inte kan samplas uniformt.
                        dx_compton, dy_compton, dz_compton = ny_steg_transformera_koordinatsystem_3d(steglängd, phi,
                                                                                                     theta,
                                                                                                     steglängd_compton,
                                                                                                     phi_compton,
                                                                                                     theta_compton)

                        # Ta ett nytt steg.
                        x, y, z, x_round, y_round, z_round = förflyttning(x, y, z, dx_compton, dy_compton, dz_compton)

                        # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                        attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size,
                                                                          z_size,
                                                                          utanför_fantom, slicad_fantom_matris,
                                                                          foton_energi)

                        if attenuerad == 0:
                            voxelvärde = slicad_fantom_matris[x_round, y_round, z_round]

                            vxv_sker = bestäm_om_vxv(voxelvärde, foton_energi, mu_max, df_attenueringsdata,
                                                     df_anatomidefinitioner)
                        else:
                            vxv_sker = True

                    # Ingångsvärden till koordinat-transformeringen (om nästa växelverkan är Comptonspridning eller Rayleighspridning).
                    theta, phi = theta_compton, phi_compton
                    steglängd = steglängd_compton

                #   ----------------------------------------------------------------------
                #   Rayleighspridning.
                #   ----------------------------------------------------------------------
                elif vxv == 'rayleigh':

                    # Sampla spridningsvinklar och steglängd.
                    theta_rayleigh = np.arccos(-1 + 2 * np.random.rand())  # ERSÄTT MED THOMSON TVÄRSNITT
                    phi_rayleigh = 2 * pi * np.random.rand()

                    vxv_sker = False
                    while vxv_sker == False:
                        # Sampla medelvägslängden från inverstransformerad attenueringsfunktion.
                        steglängd_rayleigh = medelvägslängd(mu_max)

                        # Koordinattransformation, eftersom spridningsvinklarna för Rayleighspridning inte kan samplas uniformt.
                        dx_rayleigh, dy_rayleigh, dz_rayleigh = ny_steg_transformera_koordinatsystem_3d(
                            steglängd,
                            phi,
                            theta,
                            steglängd_rayleigh,
                            phi_rayleigh,
                            theta_rayleigh)

                        # Ta ett nytt steg.
                        x, y, z, x_round, y_round, z_round = förflyttning(x, y, z, dx_rayleigh, dy_rayleigh,
                                                                          dz_rayleigh)

                        # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                        attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size,
                                                                          z_size,
                                                                          utanför_fantom, slicad_fantom_matris,
                                                                          foton_energi)

                        if attenuerad == 0:
                            voxelvärde = slicad_fantom_matris[x_round, y_round, z_round]

                            vxv_sker = bestäm_om_vxv(voxelvärde, foton_energi, mu_max, df_attenueringsdata,
                                                     df_anatomidefinitioner)
                        else:
                            vxv_sker = True

                    # Ingångsvärden till koordinat-transformeringen (om nästa växelverkan är Comptonspridning eller Rayleighspridning).
                    theta, phi = theta_rayleigh, phi_rayleigh
                    steglängd = steglängd_rayleigh

    # Några print statements, för att hålla reda på vad som händer när koden kör.
    print(f'max värdet av matrisen: {np.max(benmärg_matris_deponerad_energi)}')
    print(f'utanför: {utanför_fantom}')
    print(f'foto: {vxv_foto}')
    print(f'träffar: {träff_benmärg}')

    return benmärg_matris_deponerad_energi


def inputs_riktig_körning():
    print(
        '\n----------------------------------------------------------------------\nDags att ange parametrar.\n----------------------------------------------------------------------\n')
    print('Standard eller inte?')
    input_standard = input('\nOm standard: s, Annars: vad som helst: ')

    if input_standard == 's':
        antal_cores = 8
        iterationer_dummy = 10 ** 3
        iterationer_tot = 10 ** 5

    else:
        print(
            '\n----------------------------------------------------------------------\nVIKTIGT:\n----------------------------------------------------------------------\nAnge antal processor kärnor')
        input_antal_cores = input('Antal kärnor: ')

        if eval(input_antal_cores) > 8:
            antal_cores = 1
        else:
            antal_cores = eval(input_antal_cores)

        print(
            '\n----------------------------------------------------------------------\nDUMMY:\n----------------------------------------------------------------------\nAnge magnitud: ex 3 -> 10^3 iterationer')
        input_dummy_magnitud_iterationer = input('Magnitud: ')

        if eval(input_dummy_magnitud_iterationer) > 4:
            iterationer_dummy = 10 ** 3
        else:
            iterationer_dummy = 10 ** (eval(input_dummy_magnitud_iterationer))

        print(
            '\n----------------------------------------------------------------------\nRIKTIG:\n----------------------------------------------------------------------\nAnge skalär och magnitud: ex 5 och 5 -> 5 * 10^5 iterationer')
        input_riktig_skalär_iterationer = input('Skalär: ')
        input_riktig_magnitud_iterationer = input('Magnitud: ')

        if eval(input_riktig_magnitud_iterationer) >= 8:
            iterationer_tot = 10 ** 3
        else:
            iterationer_tot = eval(input_riktig_skalär_iterationer) * 10 ** (
                eval(input_riktig_magnitud_iterationer))

    print('antal_cores, iterationer_dummy, iterationer_tot: ', antal_cores, iterationer_dummy, iterationer_tot)
    return antal_cores, iterationer_dummy, iterationer_tot


def spara_resultat(matris, json_object):
    print(
        '\n----------------------------------------------------------------------\nDags att sparan resultaten i en matris.\n----------------------------------------------------------------------\n')
    print(
        'Om du anger namn: skriver du "text" utan citattecken kommer filerna \n[resultat_text.npy] och [inputs_text.json] \nskapas.')
    input_spara_resultat = input('\nVar vill du spara matrisen? Om standard: s, Annars: ange namn: ')

    if input_spara_resultat == 's':
        fil_namn_npy = 'resultat_multiprocess.npy'
        # Spara resultatmatrisen i en numpy fil, som sedan går att visualisera i en separat fil.
        np.save(fil_namn_npy, matris)

        fil_namn_json = 'inputs_resultat_multiprocess.json'
        with open(fil_namn_json, 'w') as f:
            f.write(json_object)
            f.close()
    else:
        fil_namn_npy = 'resultat_' + input_spara_resultat + '.npy'
        np.save(fil_namn_npy, matris)

        fil_namn_json = 'inputs_' + input_spara_resultat + '.json'
        with open(fil_namn_json, 'w') as f:
            f.write(json_object)
            f.close()

    return print(f'Resultat sparade i filerna [{fil_namn_npy}] och [{fil_namn_json}]')


if __name__ == "__main__":
    radionuklid_energi = Lu177_energi
    radionuklid_intensitet = Lu177_intensitet
    radionuklid_sannolikhet = Lu177_sannolikhet

    antal_cores, iterationer_dummy, iterationer_tot = inputs_riktig_körning()

    dictionary = {
        "antal_cores": antal_cores,
        "iterationer_dummy": iterationer_dummy,
        "iterationer_tot": iterationer_tot
    }

    json_object = json.dumps(dictionary)

    #   ----------------------------------------------------------------------
    #   Dummy run - för att snabba på den riktiga körningen av koden.
    #   ----------------------------------------------------------------------
    print(
        '\n----------------------------------------------------------------------\nDUMMY RUN\n----------------------------------------------------------------------\n')
    start = time.time()

    # Lite kod för att dela upp arbetet i flera processer, fördelat på olika processor-kärnor.
    chunk_storlek = iterationer_dummy // antal_cores
    chunk_ranges = [(i * chunk_storlek, (i + 1) * chunk_storlek) for i in range(antal_cores)]
    chunk_ranges[-1] = (chunk_ranges[-1][0], iterationer_dummy)

    # Inbakade argument för funktionen.
    args_packed = [(start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
                    slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi,
                    radionuklid_intensitet, radionuklid_sannolikhet) for start, end in chunk_ranges]

    with mp.Pool(antal_cores) as pool:
        partial_results = pool.map(run_MC_multiprocess, args_packed)

    end_time(start)

    #   ----------------------------------------------------------------------
    #   Riktig körning.
    #   ----------------------------------------------------------------------
    print(
        '\n----------------------------------------------------------------------\nACTUAL RUN\n----------------------------------------------------------------------\n')
    start = time.time()

    # Samma multiprocess-kod som ovan.
    chunk_storlek = iterationer_tot // antal_cores
    chunk_ranges = [(i * chunk_storlek, (i + 1) * chunk_storlek) for i in range(antal_cores)]
    chunk_ranges[-1] = (chunk_ranges[-1][0], iterationer_tot)

    # Samma argument, skillnaden är antalet iterationer (fotoner).
    args_packed = [(start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
                    slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi,
                    radionuklid_intensitet, radionuklid_sannolikhet) for start, end in chunk_ranges]

    with mp.Pool(antal_cores) as pool:
        partial_results = pool.map(run_MC_multiprocess, args_packed)

    # Summera resultatmatriserna från respektive process.
    benmärg_matris_deponerad_energi = np.sum(partial_results, axis=0)

    end_time(start)

    spara_resultat(benmärg_matris_deponerad_energi, json_object)

"""
----------------------------------------------------------------------
----------------------------------------------------------------------
VISUALISERING AV RESULTATEN
Visualiserar resultatmatrin, samt skapar en kapad resultatmatris.
Beräknar absorberad dos per sönderfall.
----------------------------------------------------------------------
----------------------------------------------------------------------
"""


def visualisera_resultat():
    print(
        '\n----------------------------------------------------------------------\nVilka resultat vill du visualisera?\n----------------------------------------------------------------------\n')
    print(
        'Om du anger namn: skriver du "text" utan citattecken kommer filerna \n[resultat_text.npy] och [inputs_text.json] \nvisualiseras.')
    input_resultat = input('\nOm senaste körning: s, Annars ange namn på filerna: ')

    if input_resultat == 's':
        fil_namn_npy = 'resultat_multiprocess.npy'

        resultat_matris = np.load(fil_namn_npy)
        visualisera = visualisera_matris(resultat_matris, visa_något=True)
        visualisera.show()

        slicad_resultat_matris = resultat_matris[105:160, 80:, 140:630]
        visualisera = visualisera_matris(slicad_resultat_matris, visa_något=True)
        visualisera.show()

        fil_namn_json = 'inputs_upg1_multiprocess.json'

        with (open(fil_namn_json, 'r') as f):
            json_object = json.load(f)
            iterationer_tot = json_object['iterationer_tot']
            print('\niterationer_tot: ', iterationer_tot)
            f.close()

    else:
        fil_namn_npy = 'resultat_' + input_resultat + '.npy'

        resultat_matris = np.load(fil_namn_npy)
        visualisera = visualisera_matris(resultat_matris, visa_något=True)
        visualisera.show()

        slicad_resultat_matris = resultat_matris[105:160, 80:, 140:630]
        visualisera = visualisera_matris(slicad_resultat_matris, visa_något=True)
        visualisera.show()

        fil_namn_slicad_npy = 'resultat_' + input_resultat + '_slicad.npy'

        np.save(fil_namn_slicad_npy, slicad_resultat_matris)

        fil_namn_json = 'inputs_' + input_resultat + '.json'

        with (open(fil_namn_json, 'r') as f):
            json_object = json.load(f)
            iterationer_tot = json_object['iterationer_tot']
            print('\niterationer_tot: ', iterationer_tot)
            f.close()

    return print(
        f'BENMÄRG: Energideponering per foton (eV / sönderfall): {np.sum(resultat_matris) / iterationer_tot:.1f} \nRYGGRAD: Energideponering per fotoN (eV / sönderfall): {np.sum(slicad_resultat_matris) / iterationer_tot:.1f}')


def visualisera_resultat_slicade():
    print(
        '\n----------------------------------------------------------------------\nVilka slicade resultat vill du visualisera?\n----------------------------------------------------------------------\n')
    print(
        'Om du anger namn: skriver du "text" utan citattecken kommer filerna \n[resultat_text_slicade.npy] och [inputs_text.json] \nvisualiseras.')
    input_resultat = input('\nAnge namn på filerna: ')

    fil_namn_npy = 'resultat_' + input_resultat + '_slicad.npy'

    resultat_matris = np.load(fil_namn_npy)
    visualisera = visualisera_matris(resultat_matris, visa_något=True)
    visualisera.show()

    fil_namn_json = 'inputs_' + input_resultat + '.json'

    with (open(fil_namn_json, 'r') as f):
        json_object = json.load(f)
        iterationer_tot = json_object['iterationer_tot']
        print('\niterationer_tot: ', iterationer_tot)
        f.close()

    return print(
        f'BENMÄRG: Energideponering per foton (eV / sönderfall): {np.sum(resultat_matris) / iterationer_tot:.1f} \nRYGGRAD: Energideponering per fotoN (eV / sönderfall): {np.sum(resultat_matris) / iterationer_tot:.1f}')


#   ----------------------------------------------------------------------
#   VISUALISERA RESULTATEN FRÅN upg1_sammanställt
#   ----------------------------------------------------------------------
if __name__ == "__main__":
    # visualisera_resultat()
    visualisera_resultat_slicade()
