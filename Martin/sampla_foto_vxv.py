from imports import *

@jit(nopython=True)
def foto_vxv(foton_energi):
    """
    Funktion som samplar utfallet av fotoväxelverkan.
    Antingen sker fluorescens eller inte, Augerelektroner antas deponera energi lokalt.
    :return: Energideponeringen som konsekvens av växelverkan, samt ifall fotonen slutar existera.
    """
    R = np.random.rand()

    # Fluorescens.
    if R < fluorescence_yield:
        # X-ray
        energi_deponering = foton_energi - K_alpha
        attenuerad = 0
    else:
        # Augerelektron, lokal energideponering -> finns ingen foton att följa.
        energi_deponering = foton_energi
        attenuerad = 1

    return energi_deponering, attenuerad
