from imports_file import *
from starting_position_photon import sliced_array_njure

x_size, y_size, z_size = sliced_array_njure.shape

deposited_energy = np.zeros((x_size, y_size, z_size))

