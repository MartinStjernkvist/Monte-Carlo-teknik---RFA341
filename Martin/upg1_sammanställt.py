from imports import *

from sampla_energi_start import energi_start
from matris_njure import sliced_array_njure
from sampla_position_start import position_start
from sampla_riktning_start import riktning_start
from sampla_steglängd import medelvägslängd
from sampla_steg_start import steg
from sampla_växelverkan import växelverkan


def run_MC(iterationer, mu):
    # start: sampla position, riktning och energi
    foton_energi = energi_start()
    x, y, z = position_start()
    theta, phi = riktning_start()

    while attenuerad == 0:
        # steglängd: sampla medelvägslängden från inverstransformerad attenueringsfunktion
        medelvägslängd = medelvägslängd(mu)

        #steg: gå steget till ny position i startriktning
        steg_till_ny_position = steg(theta, phi, medelvägslängd, position_start)

        växelverkan = växelverkan().bestäm_växelverkan(foton_energi)

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    return print('WORK IN PROGRESS')


print(växelverkan().bestäm_växelverkan(1000000))