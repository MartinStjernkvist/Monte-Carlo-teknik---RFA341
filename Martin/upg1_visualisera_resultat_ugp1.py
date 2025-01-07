from imports import *
from upg1_visualisera_bin_fil import visualisera_matris


def visualisera_resultat():
    print('\n----------------------------------------------------------------------\nVilka resultat vill du visualisera?\n----------------------------------------------------------------------\n')
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

        fil_namn_json = 'inputs_' + input_resultat + '.json'

        with (open(fil_namn_json, 'r') as f):
            json_object = json.load(f)
            iterationer_tot = json_object['iterationer_tot']
            print('\niterationer_tot: ', iterationer_tot)
            f.close()

    return print(f'BENMÄRG: Energideponering per foton (eV / sönderfall): {np.sum(resultat_matris) / iterationer_tot:.1f} \nRYGGRAD: Energideponering per fotoN (eV / sönderfall): {np.sum(slicad_resultat_matris) / iterationer_tot:.1f}')


#   ----------------------------------------------------------------------
#   VISUALISERA RESULTATEN FRÅN upg1_sammanställt
#   ----------------------------------------------------------------------
if __name__ == "__main__":

    visualisera_resultat()

    """
    with (open('inputs_upg1_multiprocess.json', 'r') as f):
        json_object = json.load(f)
        iterationer_tot = json_object['iterationer_tot']
        print('iterationer_tot: ', iterationer_tot)
        f.close()

    # Visualisera matrisen med registrerad energideponering.
    benmärg_matris_deponerad_energi = np.load('Martin/resultat_multiprocess.npy') #resultat_multiprocess.npy
    visualisera = visualisera_matris(benmärg_matris_deponerad_energi, visa_något=True)
    visualisera.show()

    # Visualisera en kapad version av energideponering matrisen, där endast ryggraden tas med.
    ryggrad_matris_deponerad = benmärg_matris_deponerad_energi[105:160, 80:, 140:630]
    visualisera = visualisera_matris(ryggrad_matris_deponerad, visa_något=True)
    visualisera.show()

    # Beräkna energideponering per foton: ev / foton.
    print(f'benmärg eV / decay: ', np.sum(benmärg_matris_deponerad_energi) / iterationer_tot)
    print(f'ryggrad benmärg eV / decay: ', np.sum(ryggrad_matris_deponerad) / iterationer_tot)
    """

    # np.save('resultat_5_E6', ryggrad_matris_deponerad)

    resultat_5_E6 = np.load('resultat_5_E6.npy')
    visualisera = visualisera_matris(resultat_5_E6, visa_något=True)
    visualisera.show()
    print(f'\nryggrad benmärg eV / decay: ', np.sum(resultat_5_E6) / (5 * 10 ** 6))
    """