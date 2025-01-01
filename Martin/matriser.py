from imports import *
from visualisera_bin_fil import array_phantom
from visualisera_bin_fil import visualisera_matris

#   ----------------------------------------------------------------------
#   CREATE SLICED ARRAYS
#   ----------------------------------------------------------------------

# sliced_array_phantom = array_phantom[30:230, 40:210, 570:900]

sliced_array_phantom = array_phantom[:,:,:]

x, y, z = sliced_array_phantom.shape

#   ----------------------------------------------------------------------
#   NJURE
#   ----------------------------------------------------------------------

sliced_array_njure = np.zeros((x, y, z))  # skapa tom matris för att sedan fylla i med värden

target_values = [17, 18, 19, 20, 21, 22, 23]  # värden för njure i "Anatomidefinitioner.xlsx"
mask = np.isin(sliced_array_phantom, target_values)
sliced_array_njure[mask] = sliced_array_phantom[mask]

#   ----------------------------------------------------------------------
#   BENMÄRG
#   ----------------------------------------------------------------------
"""
Samma grej som när vi gör matrisen med endast njurarna, skillnaden är att vi byter target_values till värden med benmärg
"""

# sliced_array_phantom = array_phantom[30:230, 40:210, 570:900]
# sliced_array_phantom = array_phantom[50:-50, 50:200, 600:1100] # denna matris inkluderar benmärg i rygg, samt bröst

sliced_array_benmärg = np.zeros((x, y, z))

target_values = [29]  # värdet för benmärg i "Anatomidefinitioner.xlsx"
mask = np.isin(sliced_array_phantom, target_values)
sliced_array_benmärg[mask] = sliced_array_phantom[mask]


if __name__ == "__main__":
    visualisera = visualisera_matris(sliced_array_njure)
    visualisera.show()

    visualisera = visualisera_matris(sliced_array_benmärg)
    visualisera.show()