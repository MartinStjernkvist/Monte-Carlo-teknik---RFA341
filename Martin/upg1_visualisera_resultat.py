from imports import *
from upg1_visualisera_bin_fil import visualisera_matris


def visualisera_resultat():
    """
    Visualisering av resultatmatrisen från simuleringskoden.
    """
    print(
        '\n-----------------------------------\nVilka resultat vill du visualisera?\n-----------------------------------\n')
    print(
        'Om du anger namn: skriver du "text" utan citattecken kommer filerna \n[resultat_text.npy] och [inputs_text.json] \nvisualiseras.')
    input_resultat = input('\nOm senaste körning: s, Annars ange namn på filerna: ')

    if input_resultat == 's':
        fil_namn_npy = 'resultat_multiprocess.npy'

        resultat_matris = np.load(fil_namn_npy)
        visualisera = visualisera_matris(resultat_matris, visa_något=True)
        visualisera.show()

        slicad_resultat_matris = resultat_matris[105:160, 80:, 140:630]
        # slicad_resultat_matris = resultat_matris[105:160, 80:, 240:390]
        visualisera = visualisera_matris(slicad_resultat_matris, visa_något=True)
        visualisera.show()

        fil_namn_json = 'inputs_resultat_multiprocess.json'

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
        # slicad_resultat_matris = resultat_matris[105:160, 80:, 240:390]
        visualisera = visualisera_matris(slicad_resultat_matris, visa_något=True)
        visualisera.show()

        fil_namn_slicad_npy = 'resultat_' + input_resultat + '_slicad.npy'
        # fil_namn_slicad_npy = 'resultat_v4_slicad.npy'
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
        '\n-----------------------------------\nVilka slicade resultat vill du visualisera?\n-----------------------------------\n')
    print(
        'Om du anger namn: skriver du "text" utan citattecken kommer filerna \n[resultat_text_slicade.npy] och [inputs_text.json] \nvisualiseras.')
    input_resultat = input('\nAnge namn på filerna: ')

    fil_namn_npy = 'resultat_' + input_resultat + '_slicad.npy'

    resultat_matris = np.load(fil_namn_npy)
    visualisera = visualisera_matris(resultat_matris, visa_något=2)
    visualisera.show()

    fil_namn_json = 'inputs_' + input_resultat + '.json'

    with (open(fil_namn_json, 'r') as f):
        json_object = json.load(f)
        iterationer_tot = json_object['iterationer_tot']
        print('\niterationer_tot: ', iterationer_tot)
        f.close()

    return print(
        f'BENMÄRG: Energideponering per foton (eV / sönderfall): {np.sum(resultat_matris) / iterationer_tot:.1f} \nRYGGRAD: Energideponering per fotoN (eV / sönderfall): {np.sum(resultat_matris) / iterationer_tot:.1f}')


def skapa_figurer(resultat_matris):

    x, y, z = resultat_matris.shape

    fig, ax = plt.subplots(figsize=(17.5, 5))
    plt.subplots_adjust(bottom=0.25)
    initial_slice_index = x // 2  # Start at the middle slice

    vmin, vmax = 0, 0.1 * np.max(resultat_matris)

    # Display the initial slice
    img = ax.imshow(resultat_matris[initial_slice_index, :, :], cmap='gray', vmin=vmin, vmax=vmax)
    ax.set_title(f'Slice {initial_slice_index}')
    # plt.colorbar(img, ax=ax)

    # Add a slider for navigating through slices
    ax_slider = plt.axes([0.1, 0.01, 0.75, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'Slice Index', 0, x - 1, valinit=initial_slice_index, valstep=1)


    # Update function for the slider
    def update(val):
        slice_index = int(slider.val)
        img.set_data(resultat_matris[slice_index, :, :])
        ax.set_title(f'Slice {slice_index}')
        fig.canvas.draw_idle()


    slider.on_changed(update)

    plt.tight_layout()
    plt.title('Energideponering i ryggkotorna: x-axeln.', fontsize=25)
    plt.savefig('xslice_resultat')
    plt.show()

    fig, ax = plt.subplots(figsize=(17.5, 3.5))
    plt.subplots_adjust(bottom=0.25)
    initial_slice_index = y // 2  # Start at the middle slice

    vmin, vmax = 0, 0.1 * np.max(resultat_matris)

    # Display the initial slice
    img = ax.imshow(resultat_matris[:, initial_slice_index, :], cmap='gray', vmin=vmin, vmax=vmax)
    ax.set_title(f'Slice {initial_slice_index}')
    # plt.colorbar(img, ax=ax)

    # Add a slider for navigating through slices
    ax_slider = plt.axes([0.1, 0.01, 0.75, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'Slice Index', 0, y - 1, valinit=initial_slice_index, valstep=1)


    # Update function for the slider
    def update(val):
        slice_index = int(slider.val)
        img.set_data(resultat_matris[:, slice_index, :])
        ax.set_title(f'Slice {slice_index}')
        fig.canvas.draw_idle()


    slider.on_changed(update)

    plt.tight_layout()
    plt.title('Energideponering i ryggkotorna: y-axeln.', fontsize=25)
    plt.savefig('yslice_resultat')
    plt.show()

#   -----------------------------------
#   Visualisera resultaten från simuleringskoden.
#   -----------------------------------
if __name__ == "__main__":
    visualisera_resultat()

    # visualisera_resultat_slicade()

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

    """
    # np.save('resultat_5_E6', ryggrad_matris_deponerad)

    resultat_5_E6 = np.load('resultat_5_E6.npy')
    visualisera = visualisera_matris(resultat_5_E6, visa_något=True)
    visualisera.show()
    print(f'\nryggrad benmärg eV / decay: ', np.sum(resultat_5_E6) / (5 * 10 ** 6))
    """

resultat_matris = np.load('resultat_v3_slicad.npy')
skapa_figurer(resultat_matris)
