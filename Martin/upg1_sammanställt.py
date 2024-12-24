from imports import *



# FOTON ENERGI START
from sampla_energi_start import sampla_foton_energi

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

    foton_energi = sampla_foton_energi()
    position_start = position_start()
    theta, phi = riktning_start()

    medelvägslängd = medelvägslängd(mu)
    steg = steg(theta, phi, medelvägslängd, position_start)



    return print('WORK IN PROGRESS')


