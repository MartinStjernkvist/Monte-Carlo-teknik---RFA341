from imports import *

def energi_start(radionuklid_energi, radionuklid_intensitet, radionuklid_sannolikhet):
    slump_tal = np.random.rand()
    if slump_tal <= Lu177_sannolikhet[0]:
        foton_energi = Lu177_energi[0]
    elif slump_tal <= Lu177_sannolikhet[1]:
        foton_energi = Lu177_energi[1]
    elif slump_tal <= Lu177_sannolikhet[2]:
        foton_energi = Lu177_energi[2]
    elif slump_tal <= Lu177_sannolikhet[3]:
        foton_energi = Lu177_energi[3]
    elif slump_tal <= Lu177_sannolikhet[4]:
        foton_energi = Lu177_energi[4]
    else:
        foton_energi = Lu177_energi[5]
    return foton_energi


Lu177_energi = [208_366, 112_950, 321_316, 249_674, 71_642, 136_725]  # keV

# Intensitet i % för energierna
Lu177_intensitet = [10.38, 6.2, 0.216, 0.2012, 0.1726, 0.047]

Lu177_sannolikhet = np.zeros(len(Lu177_intensitet))
for i in range(len(Lu177_intensitet)):
    if i == 0:
        Lu177_sannolikhet[i] = Lu177_intensitet[i] / np.sum(Lu177_intensitet)
    else:
        Lu177_sannolikhet[i] = Lu177_intensitet[i] / np.sum(Lu177_intensitet) + Lu177_sannolikhet[i - 1]


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
