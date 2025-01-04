from imports import *

@jit(nopython=True)
def foto_vxv(foton_energi):
    R = np.random.rand()

    if R < fluorescence_yield:
        # X-ray
        energi_deponering = foton_energi - K_alpha
        attenuerad = 0
    else:
        # augerelektron, lokal energideponering
        energi_deponering = foton_energi
        attenuerad = 1

    return energi_deponering, attenuerad
