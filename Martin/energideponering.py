from imports import *
from njure_matris import sliced_array_njure
from benmärg_matris import sliced_array_benmärg

x_size, y_size, z_size = sliced_array_njure.shape


deposited_energy = np.zeros((x_size, y_size, z_size))

