from imports import *

from sampla_energi_start import energi_start, Lu177_sannolik
from matriser import slicad_fantom_matris, slicad_njure_matris, slicad_benmärg_matris
from sampla_position_start import position_start
from sampla_riktning_och_steg_start import riktning_start, första_steg
from sampla_steglängd import medelvägslängd
from sampla_växelverkan import växelverkan
from transformation_3d import transformera_koordinatsystem
from visualisera_bin_fil import visualisera_matris
from attenueringsdata import attenueringsdata, attenueringsdata_file, anatomidefinitioner_file


def run_MC(iterationer):
    x_size, y_size, z_size = slicad_fantom_matris.shape
    benmärg_matris_deponerad_energi = np.zeros((x_size, y_size, z_size))

    utanför_fantom = 0
    vxv_foto = 0
    vxv_compton = 0
    vxv_rayleigh = 0
    träff = 0

    for i in range(iterationer):
        attenuerad = 0
        # print(i)

        # start: sampla position, riktning och energi
        foton_energi = energi_start()
        print(foton_energi)
        x_start, y_start, z_start = position_start(slicad_njure_matris)
        theta, phi = riktning_start()

        voxel_värde = slicad_fantom_matris[x_start, y_start, z_start]
        instans = attenueringsdata(voxel_värde, foton_energi, attenueringsdata_file, anatomidefinitioner_file)
        mu = instans.mu()

        # steglängd: sampla medelvägslängden från inverstransformerad attenueringsfunktion
        steglängd = medelvägslängd(mu)
        # print(f'steglängd: {steglängd}')

        # steg: gå steget till ny position i startriktning
        x, y, z = första_steg(theta, phi, steglängd, x_start, y_start, z_start)
        x_round, y_round, z_round = round(x), round(y), round(z)

        if (
                x_round < 0
                or x_round >= x_size
                or y_round < 0
                or y_round >= y_size
                or z_round < 0
                or z_round >= z_size
        ):
            utanför_fantom += 1
            i += 1

        else:
            while attenuerad == 0:

                voxel_värde = slicad_fantom_matris[x_round, y_round, z_round]
                instans = attenueringsdata(voxel_värde, foton_energi, attenueringsdata_file, anatomidefinitioner_file)
                mu = instans.mu()

                if slicad_fantom_matris[x_round, y_round, z_round] == 0:
                    utanför_fantom += 1
                    attenuerad = 1

                else:
                    instans = växelverkan(foton_energi, tvärsnitt_file)
                    vxv = instans.bestäm_växelverkan()
                    print(f'energi: {foton_energi * 10 ** (-3)} keV, vxv: {vxv}')

                    if vxv == 'foto':
                        vxv_foto += 1

                        if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                            träff += 1

                            benmärg_matris_deponerad_energi[x_round, y_round, z_round] += foton_energi

                            print(
                                f'foto: {energideponering_compton * 10 ** (-3):.2f} keV i voxel [{round(x), round(y), round(z)}]')

                        attenuerad = 1

                    elif vxv == 'compton':
                        vxv_compton += 1
                        # sampla energideponering, vinkel
                        # theta_compton = sampla_theta_compton
                        # mu_nytt energiberoende?
                        theta_compton = 1
                        phi_compton = 2 * pi * random.rand()

                        energideponering_compton = foton_energi * 0.5  # placeholder

                        if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                            träff += 1

                            benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energideponering_compton

                            print(
                                f'compton: {energideponering_compton * 10 ** (-3):.2f} keV i voxel [{round(x), round(y), round(z)}]')

                        foton_energi = foton_energi - energideponering_compton




                        steglängd_compton = medelvägslängd(mu)
                        vektor_compton, _, dx_compton, dy_compton, dz_compton = transformera_koordinatsystem(steglängd,
                                                                                                             phi,
                                                                                                             theta,
                                                                                                             steglängd_compton,
                                                                                                             phi_compton,
                                                                                                             theta_compton)

                        x_round = round(x + dx_compton / voxel_sidlängd)
                        y_round = round(y + dy_compton / voxel_sidlängd)
                        z_round = round(z + dz_compton / voxel_sidlängd)

                        if (
                                x_round < 0
                                or x_round >= x_size
                                or y_round < 0
                                or y_round >= y_size
                                or z_round < 0
                                or z_round >= z_size
                        ):
                            utanför_fantom += 1
                            attenuerad = 1

                        theta, phi = theta_compton, phi_compton
                        steglängd = steglängd_compton

                    elif vxv == 'rayleigh':
                        vxv_rayleigh += 1
                        theta_rayleigh = 1
                        phi_rayleigh = 2 * pi * random.rand()
                        steglängd_rayleigh = medelvägslängd(mu)

                        vektor_rayleigh, _, dx_rayleigh, dy_rayleigh, dz_rayleigh = transformera_koordinatsystem(
                            steglängd,
                            phi,
                            theta,
                            steglängd_rayleigh,
                            phi_rayleigh,
                            theta_rayleigh)

                        x_round = round(x + dx_rayleigh / voxel_sidlängd)
                        y_round = round(y + dy_rayleigh / voxel_sidlängd)
                        z_round = round(z + dz_rayleigh / voxel_sidlängd)

                        if (
                                x_round < 0
                                or x_round >= x_size
                                or y_round < 0
                                or y_round >= y_size
                                or z_round < 0
                                or z_round >= z_size
                        ):
                            utanför_fantom += 1
                            attenuerad = 1

                        theta, phi = theta_rayleigh, phi_rayleigh
                        steglängd = steglängd_rayleigh

    print(f'max värdet av matrisen: {np.max(benmärg_matris_deponerad_energi)}')
    print(f'utanför: {utanför_fantom}')
    print(f'foto: {vxv_foto}')
    print(f'rayleigh: {vxv_rayleigh}')
    print(f'compton: {vxv_compton}')
    print(f'träffar: {träff}')
    return benmärg_matris_deponerad_energi


if __name__ == "__main__":
    start = time.time()

    tvärsnitt_file = '../given_data/Tvärsnittstabeller_Fotoner.xlsx'
    df = pd.read_excel(tvärsnitt_file, index_col=None)

    iterationer = 100
    benmärg_matris_deponerad_energi = run_MC(iterationer)

    visualisera = visualisera_matris(benmärg_matris_deponerad_energi)
    end_time(start)

    visualisera.show()

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
