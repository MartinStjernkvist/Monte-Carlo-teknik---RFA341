from imports import *


@jit(nopython=True)
def compton_vinkel_och_energiförlust(foton_energi):
    """
    Khans rejektionsalgoritm, artikeln "A Monte Carlo program for the simulation of scintillation camera characteristics"
    Algoritmen samplar vinkeln och energideponeringen för en foton som Comptonväxelverkar.
    :return: Spridningsvinkeln theta, den spridda fotonens energi och den spridda fotonens energiförlust (=energideponering).
    """

    while True:
        # Slumpa tre tal.
        R_1 = np.random.rand()
        R_2 = np.random.rand()
        R_3 = np.random.rand()

        # Termen alpha är fotonens energi i termer av vilomassan (energi) för en elektron.
        alpha = foton_energi / E_e

        if R_1 <= (2 * alpha + 1) / (2 * alpha + 9):

            eta = 2 * alpha * R_2
            if R_3 <= 4 * (eta ** -1 - eta ** -2):
                theta = np.arccos(1 - 2 * R_2)
                spridd_foton_energi = foton_energi / eta
                energi_förlust = foton_energi - spridd_foton_energi

                break

        elif R_1 > (2 * alpha + 1) / (2 * alpha + 9):

            eta = (2 * alpha + 1) / (2 * R_2 * alpha + 1)
            theta = np.arccos(1 - (eta - 1) / alpha)

            if R_3 <= 0.5 * (np.cos(theta) ** 2 + eta ** -1):
                spridd_foton_energi = foton_energi / eta
                theta = np.arccos(1 - (eta - 1) / alpha)
                energi_förlust = foton_energi - spridd_foton_energi

                break

        else:
            continue

    return theta, spridd_foton_energi, energi_förlust
