from imports import *

from sampla_energi_start import energi_start
from matriser import sliced_array_phantom, slicad_njure_matris, slicad_benmärg_matris
from sampla_position_start import position_start
from sampla_riktning_start import riktning_start
from sampla_steglängd import medelvägslängd
from sampla_steg_start import första_steg
from sampla_växelverkan import växelverkan
from transformation_3d import transformera_koordinatsystem
from energideponering import benmärg_matris_deponerad_energi, energideponering_benmärg


def run_MC(iterationer, mu):
    x_size, y_size, z_size = sliced_array_phantom.shape

    # start: sampla position, riktning och energi
    foton_energi = energi_start()
    x_start, y_start, z_start = position_start()
    theta, phi = riktning_start()

    tot_x = x_start
    tot_y = y_start
    tot_z = z_start

    while attenuerad == 0:
        # steglängd: sampla medelvägslängden från inverstransformerad attenueringsfunktion
        steglängd = medelvägslängd(mu)

        # steg: gå steget till ny position i startriktning
        x, y, z = första_steg(theta, phi, steglängd, x_start, y_start, z_start)

        if x >= x_size or y >= y_size or z >= z_size or sliced_array_phantom[x, y, z] == 0:
            attenuerad = 1

        if slicad_benmärg_matris[x, y, z] == 0:
            continue

        växelverkan = växelverkan().bestäm_växelverkan(foton_energi)

        if växelverkan == 'foto':
            attenuerad = 1
            energideponering_benmärg(x, y, z, foton_energi)

        if växelverkan == 'compton':
            # sampla energideponering, vinkel
            # theta_compton = sampla_theta_compton
            # mu_nytt energiberoende?
            theta_compton = 1
            phi_compton = 2 * pi * random.rand()
            mu_ny = 1
            energideponering_compton = 1

            energideponering_benmärg(x, y, z, energideponering_compton)

            steglängd_compton = medelvägslängd(mu_ny)
            vektor_compton, _ = transformera_koordinatsystem(steglängd, phi, theta, steglängd_compton, phi_compton,
                                                             theta_compton)

            tot_x += vektor_compton[0]
            tot_y += vektor_compton[1]
            tot_z += vektor_compton[2]

            x_start = x + vektor_compton[0]
            y_start = y + vektor_compton[1]
            z_start = z + vektor_compton[2]

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    return benmärg_matris_deponerad_energi


print(växelverkan().bestäm_växelverkan(1000000))
