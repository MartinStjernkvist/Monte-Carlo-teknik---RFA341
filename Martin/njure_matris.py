from imports import *
from visualisera_bin_fil import array_phantom
from visualisera_bin_fil import visualisera_matris

#   ----------------------------------------------------------------------
#   CREATE SLICED ARRAYS
#   ----------------------------------------------------------------------

sliced_array_phantom = array_phantom[30:230, 40:210, 570:900]
# sliced_array_phantom = array_3d[:,:,:]
x, y, z = sliced_array_phantom.shape

sliced_array_njure = np.zeros((x, y, z)) # skapa tom matris för att sedan fylla i med värden

target_values = [17, 18, 19, 20, 21, 22, 23] # värden för njure i "Anatomidefinitioner.xlsx"
mask = np.isin(sliced_array_phantom, target_values)
sliced_array_njure[mask] = sliced_array_phantom[mask]

visualisera = visualisera_matris(sliced_array_njure)
visualisera.show()