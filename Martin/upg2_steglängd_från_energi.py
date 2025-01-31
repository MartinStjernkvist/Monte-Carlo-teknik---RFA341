from imports import *
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd


def steglängd_från_energi(energi_start, rho_medium, stopping_power_data, energiförlust_faktor):
    """
    Funktion som beräknar steglängden utifrån startenergi och energiförlust.
    :param energi_start: Startenergin för elektronen.
    :param rho_medium: Mediumets densitet (kg/m^3).
    :param stopping_power_data: Stopping power tabelldata.
    :param energiförlust_faktor: Faktor för energiförlust, t.ex. 0.97 för 3%.
    :return: Steglängden till nästa kollision.
    """

    energi_1 = energi_start

    # Bestäm ny energi, varefter steglängden kan lösas ut m.h.a. stopping power.
    energi_2 = energi_1 * energiförlust_faktor  # antag energiförlust innan nästa kollission

    # Stopping power från tabelldata, med hjälp av linjärinterpolering.
    # Delvis stopping power för energi_1, delvis för energi_2.
    stopping_power_1, _ = stopping_power_och_steglängd(energi_1, rho_medium, stopping_power_data)
    stopping_power_2, _ = stopping_power_och_steglängd(energi_2, rho_medium, stopping_power_data)

    # Medelvärdet av stopping power mellan startpunkt och slutpunkt.
    stopping_power_medel = (stopping_power_2 + stopping_power_1) / 2
    # print(f'stopping power medel: {stopping_power_medel * 10**(-6):.3f} MeV / m')

    # Beräkna steglängden.
    steglängd = - (energi_2 - energi_1) / stopping_power_medel
    return steglängd


if __name__ == "__main__":

    rho_medium = rho_vatten
    energiförlust_faktor = 0.97
    def steglängd_från_energi_kort(energi_start):
        """
        Funktion som beräknar steglängden utifrån startenergi och energiförlust.
        :param energi_start: Startenergin för elektronen.
        :param rho_medium: Mediumets densitet (kg/m^3).
        :param stopping_power_data: Stopping power tabelldata.
        :param energiförlust_faktor: Faktor för energiförlust, t.ex. 0.97 för 3%.
        :return: Steglängden till nästa kollision.
        """

        energi_1 = energi_start

        # Bestäm ny energi, varefter steglängden kan lösas ut m.h.a. stopping power.
        energi_2 = energi_1 * energiförlust_faktor  # antag energiförlust innan nästa kollission

        # Stopping power från tabelldata, med hjälp av linjärinterpolering.
        # Delvis stopping power för energi_1, delvis för energi_2.
        stopping_power_1, _ = stopping_power_och_steglängd(energi_1, rho_medium, stopping_power_data)
        stopping_power_2, _ = stopping_power_och_steglängd(energi_2, rho_medium, stopping_power_data)

        # Medelvärdet av stopping power mellan startpunkt och slutpunkt.
        stopping_power_medel = (stopping_power_2 + stopping_power_1) / 2
        # print(f'stopping power medel: {stopping_power_medel * 10**(-6):.3f} MeV / m')

        # Beräkna steglängden.
        steglängd = - (energi_2 - energi_1) / stopping_power_medel
        return steglängd


    stopping_power_data = np.loadtxt(elektron_stopping_power_data)

    energi_list = np.linspace(50, 2_280_000)
    steglängd_list = list(map(steglängd_från_energi_kort, energi_list))
    plt.plot(energi_list, steglängd_list)

    # Visa figuren
    plt.show()