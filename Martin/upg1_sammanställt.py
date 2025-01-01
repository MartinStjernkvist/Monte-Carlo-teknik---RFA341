from imports import *

from sampla_energi_start import energi_start
from matriser import sliced_array_phantom, sliced_array_njure, sliced_array_benmärg
from sampla_position_start import position_start
from sampla_riktning_start import riktning_start
from sampla_steglängd import medelvägslängd
from sampla_steg_start import första_steg
from sampla_växelverkan import växelverkan
from transformation_3d import transformera_koordinatsystem


def run_MC(iterationer, mu):
    x_size, y_size, z_size = sliced_array_phantom.shape

    # start: sampla position, riktning och energi
    foton_energi = energi_start()
    x_start, y_start, z_start = position_start()
    theta, phi = riktning_start()


    while attenuerad == 0:
        # steglängd: sampla medelvägslängden från inverstransformerad attenueringsfunktion
        medelvägslängd = medelvägslängd(mu)

        #steg: gå steget till ny position i startriktning
        x,y,z = första_steg(theta, phi, medelvägslängd, x_start, y_start, z_start)

        if x>= x_size or y >= y_size or z >= z_size or sliced_array_phantom[x,y,z] == 0:
            attenuerad = 1

        if sliced_array_benmärg[x,y,z] == 0:
            continue

        växelverkan = växelverkan().bestäm_växelverkan(foton_energi)

        if växelverkan == 'foto':
            attenuerad = 1


        if växelverkan == 'compton':
            # sampla energideponering, vinkel
            theta = 1
            phi = 1
            x_start = x
            y_start = y
            z_start = z



    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    return print('WORK IN PROGRESS')


print(växelverkan().bestäm_växelverkan(1000000))