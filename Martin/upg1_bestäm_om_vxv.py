from imports import *
from upg1_attenueringsdata import attenueringsdata

def bestäm_om_vxv(voxelvärde, energi, mu_max, df_attenueringsdata, df_anatomidefinitioner):

    instans = attenueringsdata(voxelvärde, energi, df_attenueringsdata, df_anatomidefinitioner)
    mu = instans.mu()

    if np.random.rand() <= mu / mu_max:
        vxv_sker = True
    else:
        vxv_sker = False

    return vxv_sker