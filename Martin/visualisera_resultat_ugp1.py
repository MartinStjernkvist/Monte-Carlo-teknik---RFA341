from imports import *
from visualisera_bin_fil import visualisera_matris
from input_upg1_multiprocess import iterationer_tot

#   ----------------------------------------------------------------------
#   VISUALISERA RESULTATEN FRÅN upg1_sammanställt
#   ----------------------------------------------------------------------

if __name__ == "__main__":

    # Visualisera matrisen med registrerad energideponering.
    benmärg_matris_deponerad_energi = np.load('resultat_multiprocess.npy')
    visualisera = visualisera_matris(benmärg_matris_deponerad_energi, visa_något=True)
    visualisera.show()

    # Visualisera en kapad version av energideponering matrisen, där endast ryggraden tas med.
    ryggrad_matris_deponerad = benmärg_matris_deponerad_energi[105:155, 75:, 260:630]
    visualisera = visualisera_matris(ryggrad_matris_deponerad, visa_något=True)
    visualisera.show()

    # Beräkna energideponering per foton: ev / foton.
    print(f'benmärg eV / decay: ', np.sum(benmärg_matris_deponerad_energi) / iterationer_tot)
    print(f'ryggrad benmärg eV / decay: ', np.sum(ryggrad_matris_deponerad) / iterationer_tot)