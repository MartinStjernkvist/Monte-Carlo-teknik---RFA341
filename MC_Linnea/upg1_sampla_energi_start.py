from imports import *


def energi_start(radionuklid_energi, radionuklid_sannolikhet):
    """
    Funktion som samplar den ursprungliga partikelenergin för varje ny partikel som sänds ut.
    :param radionuklid_energi: Lista med partikelenergier som är möjliga.
    :param radionuklid_sannolikhet: Lista med sannolikheten för partikelergierna.
    :return:
    """
    slump_tal = np.random.rand()
   
    if slump_tal <= radionuklid_sannolikhet[0]:
        partikel_energi = radionuklid_energi[0]
    elif slump_tal <= radionuklid_sannolikhet[1]:
        partikel_energi = radionuklid_energi[1]
    elif slump_tal <= radionuklid_sannolikhet[2]:
        partikel_energi = radionuklid_energi[2]
    elif slump_tal <= radionuklid_sannolikhet[3]:
        partikel_energi = radionuklid_energi[3]
    elif slump_tal <= radionuklid_sannolikhet[4]:
        partikel_energi = radionuklid_energi[4]
    else:
        partikel_energi = radionuklid_energi[5]
    """
    for i in range(len(radionuklid_sannolikhet)):
        if slump_tal <= radionuklid_sannolikhet[i]:
            partikel_energi = radionuklid_energi[i]
        else:
            i+=1
    """
    return partikel_energi


