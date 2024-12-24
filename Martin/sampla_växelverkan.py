from imports import *

def sigma_compton(energi):
    return print('WIP')

def sigma_foto(energi):
    return print('WIP')

def sigma_par(energi):
    return print('WIP')

def växelverkan(energi):
    sigma_compton = sigma_compton(energi)
    sigma_foto = sigma_foto(energi)
    sigma_par = sigma_par(energi)

    Lu177_energ = [sigma_foto, sigma_compton, sigma_par]  # keV

    # Intensitet i % för energierna
    Lu177_intens = [10.38, 6.2, 0.216, 0.2012, 0.1726, 0.047]

    Lu177_sannolik = np.zeros(len(Lu177_intens))
    for i in range(len(Lu177_intens)):
        if i == 0:
            Lu177_sannolik[i] = Lu177_intens[i] / np.sum(Lu177_intens)
        else:
            Lu177_sannolik[i] = Lu177_intens[i] / np.sum(Lu177_intens) + Lu177_sannolik[i - 1]

    rand = np.random.rand()
