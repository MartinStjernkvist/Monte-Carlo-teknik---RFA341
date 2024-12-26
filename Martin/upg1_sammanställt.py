from imports import *



# FOTON ENERGI START
from sampla_energi_start import energi_start

# MATRIS NJURE, POSITION START
from njure_matris import sliced_array_njure
from sampla_position_start import position_start

# RIKTNING
from sampla_riktning_start import riktning_start

# MEDELVÄGSLÄNGD
from sampla_steglängd import medelvägslängd

# STEG TILL NY POSITION
from steg_start import steg

def run_MC(iterationer, mu):

    foton_energi = energi_start()
    x,y,z = position_start()
    theta, phi = riktning_start()

    medelvägslängd = medelvägslängd(mu)
    steg_till_ny_position = steg(theta, phi, medelvägslängd, position_start)

    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)


    return print('WORK IN PROGRESS')


