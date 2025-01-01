from imports import *

from sampla_energi_start import energi_start, Lu177_sannolik
from matriser import sliced_array_phantom, slicad_njure_matris, slicad_benmärg_matris
from sampla_position_start import position_start
from sampla_riktning_start import riktning_start
from sampla_steglängd import medelvägslängd
from sampla_steg_start import första_steg
from sampla_växelverkan import växelverkan
from transformation_3d import transformera_koordinatsystem
from energideponering import benmärg_matris_deponerad_energi, energideponering_benmärg
from visualisera_bin_fil import visualisera_matris


def run_MC(iterationer, mu):
    for i in range(iterationer):
        attenuerad = 0

        x_size, y_size, z_size = sliced_array_phantom.shape

        # start: sampla position, riktning och energi
        foton_energi = energi_start()
        x_start, y_start, z_start = position_start(slicad_njure_matris)
        theta, phi = riktning_start()

        # x_tot = x_start
        # y_tot = y_start
        # z_tot = z_start

        while attenuerad == 0:
            # steglängd: sampla medelvägslängden från inverstransformerad attenueringsfunktion
            steglängd = medelvägslängd(mu)

            # steg: gå steget till ny position i startriktning
            x, y, z = första_steg(theta, phi, steglängd, x_start, y_start, z_start)

            if x >= x_size or y >= y_size or z >= z_size or sliced_array_phantom[round(x), round(y), round(z)] == 0:
                attenuerad = 1

            instans = växelverkan(foton_energi, tvärsnitt_file)
            vxv = instans.bestäm_växelverkan()

            if vxv == 'foto':
                attenuerad = 1
                energideponering_benmärg(x, y, z, foton_energi, slicad_benmärg_matris)
                print(f'foto energideponering: energi {energideponering_compton} i voxel [{x:.2f}, {y:.2f}, {z:.2f}]')

            elif vxv == 'compton':
                # sampla energideponering, vinkel
                # theta_compton = sampla_theta_compton
                # mu_nytt energiberoende?
                theta_compton = 1
                phi_compton = 2 * pi * random.rand()
                mu_ny = 1
                energideponering_compton = foton_energi * 0.5  # placeholder

                energideponering_benmärg(x, y, z, energideponering_compton, benmärg_matris_deponerad_energi)
                print(
                    f'compton energideponering: energi {energideponering_compton * 10 ** (-3):.2f} keV i voxel [{x:.2f}, {y:.2f}, {z:.2f}]')
                foton_energi = foton_energi - energideponering_compton

                steglängd_compton = medelvägslängd(mu_ny)
                vektor_compton, _, dx_compton, dy_compton, dz_compton = transformera_koordinatsystem(steglängd, phi,
                                                                                                     theta,
                                                                                                     steglängd_compton,
                                                                                                     phi_compton,
                                                                                                     theta_compton)

                # x_tot += vektor_compton[0]
                # y_tot += vektor_compton[1]
                # z_tot += vektor_compton[2]

                x_start = x + dx_compton
                y_start = y + dy_compton
                z_start = z + dz_compton

                theta, phi = theta_compton, phi_compton
                steglängd = steglängd_compton

            else:
                theta_rayleigh = 1
                phi_rayleigh = 2 * pi * random.rand()
                steglängd_rayleigh = medelvägslängd(mu)

                vektor_rayleigh, _, dx_rayleigh, dy_rayleigh, dz_rayleigh = transformera_koordinatsystem(steglängd, phi,
                                                                                                         theta,
                                                                                                         steglängd_rayleigh,
                                                                                                         phi_rayleigh,
                                                                                                         theta_rayleigh)

                # x_tot += vektor_rayleigh[0]
                # y_tot += vektor_rayleigh[1]
                # z_tot += vektor_rayleigh[2]

                x_start = x + dx_rayleigh
                y_start = y + dy_rayleigh
                z_start = z + dz_rayleigh

                theta, phi = theta_rayleigh, phi_rayleigh
                steglängd = steglängd_rayleigh

        # r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    return benmärg_matris_deponerad_energi


if __name__ == "__main__":
    start = time.time()

    tvärsnitt_file = '../given_data/Tvärsnittstabeller_Fotoner.xlsx'
    df = pd.read_excel(tvärsnitt_file, index_col=None)

    iterationer = 100
    mu = 0.5
    benmärg_matris_deponerad_energi = run_MC(iterationer, mu)

    visualisera = visualisera_matris(benmärg_matris_deponerad_energi)
    end_time(start)
    
    visualisera.show()

