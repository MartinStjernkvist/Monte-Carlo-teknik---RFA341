from imports import *
from visualisera_bin_fil import visualisera_matris
from input_upg1_multiprocess import iterationer_tot

#   ----------------------------------------------------------------------
#   VISUALISERA RESULTATEN FRÅN upg1_sammanställt
#   ----------------------------------------------------------------------

if __name__ == "__main__":
    # benmärg_matris_deponerad_energi = np.load('benmärg_matris_deponerad_energi.npy')
    # visualisera = visualisera_matris(benmärg_matris_deponerad_energi, visa_något=True)
    #
    # visualisera.show()

    benmärg_matris_deponerad_energi = np.load('resultat_multiprocess.npy')
    visualisera = visualisera_matris(benmärg_matris_deponerad_energi, visa_något=True)
    visualisera.show()

    ryggrad_matris_deponerad = benmärg_matris_deponerad_energi[105:155, 75:, 260:630]
    visualisera = visualisera_matris(ryggrad_matris_deponerad, visa_något=True)
    visualisera.show()

    print(f'benmärg eV / decay: ', np.sum(benmärg_matris_deponerad_energi) / iterationer_tot)
    print(f'ryggrad benmärg eV / decay: ', np.sum(ryggrad_matris_deponerad) / iterationer_tot)

    # fig, ax = plt.subplots(figsize=(8, 8))
    # plt.subplots_adjust(bottom=0.25)
    # x_size, y_size, z_size = benmärg_matris_deponerad_energi.shape
    # initial_slice_index = z_size // 2  # Start at the middle slice
    #
    # # Display the initial slice
    # img = ax.imshow(benmärg_matris_deponerad_energi[:, :, initial_slice_index], cmap='hot')
    # ax.set_title(f'Slice {initial_slice_index}')
    # plt.colorbar(img, ax=ax)
    #
    # # Add a slider for navigating through slices
    # ax_slider = plt.axes([0.25, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    # slider = Slider(ax_slider, 'Slice Index', 0, z_size - 1, valinit=initial_slice_index, valstep=1)
    #
    #
    # # Update function for the slider
    # def update(val):
    #     slice_index = int(slider.val)
    #     img.set_data(benmärg_matris_deponerad_energi[:, :, slice_index])
    #     ax.set_title(f'Slice {slice_index}')
    #     fig.canvas.draw_idle()
    #
    #
    # slider.on_changed(update)
    #
    # plt.show()
