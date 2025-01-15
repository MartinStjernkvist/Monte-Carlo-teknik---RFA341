from imports import *
from upg1_attenueringsdata import attenueringsdata

def bestäm_om_vxv_sker(voxelvärde, energi, mu_max, df_attenueringsdata, df_anatomidefinitioner):
    """
    Funktion som bestämmer om växelverkan sker efter
    en steglängd har tagits.
    :param voxelvärde:
    :param energi:
    :param mu_max:
    :param df_attenueringsdata:
    :param df_anatomidefinitioner:
    :return:
    """

    # Attenueringskoefficienten i den nya positionen.
    instans = attenueringsdata(voxelvärde, energi, df_attenueringsdata, df_anatomidefinitioner)
    mu = instans.mu()

    if np.random.rand() <= mu / mu_max:
        vxv_sker = True
    else:
        vxv_sker = False

    return vxv_sker