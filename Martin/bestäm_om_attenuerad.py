from imports import *

@jit(nopython=True)
def bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size, utanför_fantom, slicad_fantom_matris, foton_energi):
    """
    Funktion som bestämmer om en foton ska fortsätta följas eller inte.
    Termen "attenuerad" i koden kan även vara att fotonen rymt från fantommatrisen.
    :param x_round: voxelposition x
    :param y_round: voxelposition y
    :param z_round: voxelposition z
    :param x_size: storlek matris x
    :param y_size: storlek matris y
    :param z_size: storlek matris z
    :param utanför_fantom: värde som ökar med 1, för att hålla reda på vad som händer när koden körs
    :param slicad_fantom_matris: fantommatrisen
    :return: om attenuerad = 1 kommer fotonen att sluta följas
    """

    # Om positionen för fotonen är utanför fantommatrisen, indikera att fotonen ska sluta följas.
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

    # Om voxelvärdet för positionen är innanför matrisen, men befinner sig i luft (voxelvärde = 0)
    # eller
    # fotonenergin är mindre än gränsenergin för fotoabsorption -> indikera att fotonen ska sluta följas.
    elif slicad_fantom_matris[x_round, y_round, z_round] == 0 or foton_energi < foton_energi_threshhold:
        utanför_fantom += 1
        attenuerad = 1

    else:
        attenuerad = 0

    return attenuerad, utanför_fantom
