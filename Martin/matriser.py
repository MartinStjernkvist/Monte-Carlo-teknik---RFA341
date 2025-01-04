from imports import *
from visualisera_bin_fil import fantom_matris, visualisera_matris

#   ----------------------------------------------------------------------
#   CREATE SLICED ARRAYS
#   ----------------------------------------------------------------------

# slicad_fantom_matris = fantom_matris[30:230, 40:210, 570:900]
# slicad_fantom_matris = fantom_matris[50:-50, 50:200, 600:1100] # denna matris inkluderar benmärg i rygg, samt bröst
# slicad_fantom_matris = fantom_matris
# slicad_fantom_matris = fantom_matris[50:-50, 50:200, 600:1100]  # denna matris inkluderar benmärg i rygg, samt bröst

# slicad_fantom_matris = fantom_matris[50:-50, 50:200, 650:1050]  # samma som ovan, fast har endast med intervallet varvid 10**5 iterationer ger något utslag

slicad_fantom_matris = fantom_matris[:, 34:212, 450:1182] # inkluderar fantom, avskuret vid benen

x, y, z = slicad_fantom_matris.shape

#   ----------------------------------------------------------------------
#   NJURE
#   ----------------------------------------------------------------------

slicad_njure_matris = np.zeros((x, y, z))  # skapa tom matris för att sedan fylla i med värden

target_values = [17, 18, 19, 20, 21, 22, 23]  # värden för njure i "Anatomidefinitioner.xlsx"
mask = np.isin(slicad_fantom_matris, target_values)
slicad_njure_matris[mask] = slicad_fantom_matris[mask]

#   ----------------------------------------------------------------------
#   BENMÄRG
#   ----------------------------------------------------------------------
"""
Samma grej som för njurarna, skillnaden är att vi byter target_values till värden med benmärg
"""

slicad_benmärg_matris = np.zeros((x, y, z))

target_values = [29]  # värdet för benmärg i "Anatomidefinitioner.xlsx"
mask = np.isin(slicad_fantom_matris, target_values)
slicad_benmärg_matris[mask] = slicad_fantom_matris[mask]

if __name__ == "__main__":
    visualisera = visualisera_matris(slicad_fantom_matris, visa_något=True)
    # visualisera = visualisera_matris(slicad_fantom_matris)
    visualisera.show()

    visualisera = visualisera_matris(slicad_njure_matris)
    visualisera.show()

    visualisera = visualisera_matris(slicad_benmärg_matris)
    visualisera.show()