from imports import *
from visualisera_bin_fil import array_phantom
from visualisera_bin_fil import visualisera_matris

"""
Samma grej som när vi gör matrisen med endast njurarna, skillnaden är att vi byter target_values till värden med benmärg
"""

#   ----------------------------------------------------------------------
#   CREATE SLICED ARRAYS
#   ----------------------------------------------------------------------

# sliced_array_phantom = array_phantom[30:230, 40:210, 570:900]
sliced_array_phantom = array_phantom[50:-50,50:200,600:1100] # denna matris inkluderar benmärg i rygg, samt bröst

x, y, z = sliced_array_phantom.shape

sliced_array_benmärg = np.zeros((x, y, z))

target_values = [29] # värdet för benmärg i "Anatomidefinitioner.xlsx"
mask = np.isin(sliced_array_phantom, target_values)
sliced_array_benmärg[mask] = sliced_array_phantom[mask]

visualisera = visualisera_matris(sliced_array_benmärg)
visualisera.show()