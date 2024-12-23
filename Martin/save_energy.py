from imports_file import *
from starting_position_photon import sliced_array_njure
from end_position import sliced_array_benmärg

x_size, y_size, z_size = sliced_array_njure.shape


deposited_energy = np.zeros((x_size, y_size, z_size))

