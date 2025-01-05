from imports import *


def energi_start(radionuklid_energi, radionuklid_sannolikhet):
    """
    Funktion som samplar den ursprungliga fotonenergin för varje ny foton som sänds ut.
    :param radionuklid_energi: Lista med fotonenergier som är möjliga.
    :param radionuklid_sannolikhet: Lista med sannolikheten för fotonergierna.
    :return:
    """

    slump_tal = np.random.rand()

    if slump_tal <= radionuklid_sannolikhet[0]:
        foton_energi = radionuklid_energi[0]
    elif slump_tal <= radionuklid_sannolikhet[1]:
        foton_energi = radionuklid_energi[1]
    elif slump_tal <= radionuklid_sannolikhet[2]:
        foton_energi = radionuklid_energi[2]
    elif slump_tal <= radionuklid_sannolikhet[3]:
        foton_energi = radionuklid_energi[3]
    elif slump_tal <= radionuklid_sannolikhet[4]:
        foton_energi = radionuklid_energi[4]
    else:
        foton_energi = radionuklid_energi[5]

    return foton_energi

if __name__ == "__main__":
    Antal_iterationer = 100

    Foton_energi = []

    for i in range(Antal_iterationer):
        Slump_tal = np.random.rand()
        if Slump_tal <= Lu177_sannolikhet[0]:
            Foton_energi.append(Lu177_energi[0])
        elif Slump_tal <= Lu177_sannolikhet[1]:
            Foton_energi.append(Lu177_energi[1])
        elif Slump_tal <= Lu177_sannolikhet[2]:
            Foton_energi.append(Lu177_energi[2])
        elif Slump_tal <= Lu177_sannolikhet[3]:
            Foton_energi.append(Lu177_energi[3])
        elif Slump_tal <= Lu177_sannolikhet[4]:
            Foton_energi.append(Lu177_energi[4])
        else:
            Foton_energi.append(Lu177_energi[5])

    fig = plt.figure()
    plt.plot([1, 2, 3, 4, 5, 6], Lu177_sannolikhet)
    plt.show()

    print(np.sum(Lu177_sannolikhet))
