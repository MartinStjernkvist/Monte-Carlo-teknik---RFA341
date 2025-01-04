from imports import *

@jit(nopython=True)
def bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size, utanför_fantom, slicad_fantom_matris, foton_energi):
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

    elif slicad_fantom_matris[x_round, y_round, z_round] == 0 or foton_energi < foton_energi_threshhold:
        utanför_fantom += 1
        attenuerad = 1

    else:
        attenuerad = 0

    return attenuerad, utanför_fantom
