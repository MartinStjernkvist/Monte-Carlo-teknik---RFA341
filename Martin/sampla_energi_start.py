from imports import *

if __name__ == "__main__":
    # Ursprungsenergier från nukliden
    Lu177_energ = [208.3, 112.9, 321.3, 249.7, 71.6, 136.7]  # keV

    # Intensitet i % för energierna
    Lu177_intens = [10.38, 6.2, 0.216, 0.2012, 0.1726, 0.047]

    Lu177_sannolik = np.zeros(len(Lu177_intens))
    for i in range(len(Lu177_intens)):
        if i == 0:
            Lu177_sannolik[i] = Lu177_intens[i] / np.sum(Lu177_intens)
        else:
            Lu177_sannolik[i] = Lu177_intens[i] / np.sum(Lu177_intens) + Lu177_sannolik[i - 1]

    Antal_iterationer = 100

    Foton_energi = []

    for i in range(Antal_iterationer):
        Slump_tal = np.random.rand()
        if Slump_tal <= Lu177_sannolik[0]:
            Foton_energi.append(Lu177_energ[0])
        elif Slump_tal <= Lu177_sannolik[1]:
            Foton_energi.append(Lu177_energ[1])
        elif Slump_tal <= Lu177_sannolik[2]:
            Foton_energi.append(Lu177_energ[2])
        elif Slump_tal <= Lu177_sannolik[3]:
            Foton_energi.append(Lu177_energ[3])
        elif Slump_tal <= Lu177_sannolik[4]:
            Foton_energi.append(Lu177_energ[4])
        else:
            Foton_energi.append(Lu177_energ[5])

    fig = plt.figure()
    plt.plot([1, 2, 3, 4, 5, 6], Lu177_sannolik)
    plt.show()

    print(np.sum(Lu177_sannolik))


def energi_start():
    slump_tal = np.random.rand()
    if slump_tal <= Lu177_sannolik[0]:
        foton_energi = Lu177_energ[0]
    elif slump_tal <= Lu177_sannolik[1]:
        foton_energi = Lu177_energ[1]
    elif slump_tal <= Lu177_sannolik[2]:
        foton_energi = Lu177_energ[2]
    elif slump_tal <= Lu177_sannolik[3]:
        foton_energi = Lu177_energ[3]
    elif slump_tal <= Lu177_sannolik[4]:
        foton_energi = Lu177_energ[4]
    else:
        foton_energi = Lu177_energ[5]
    return foton_energi
