from imports import *


def stopping_power_och_steglängd(energi, rho_medium, stopping_power_data):
    """
    Funktion som ger stopping power och
    :param energi: Partikelns energi i eV.
    :return: Stopping power (eV/m) och steglängden (m).
    """

    energi_MeV_list = stopping_power_data[:, 0]
    energi_list = np.array([(lambda x: x * 10**6)(x) for x in energi_MeV_list]) #Fråm MeV till eV

    stopping_power_MeV_list = stopping_power_data[:, 1]
    STP_list =np.array([(lambda x: x * 10**5)(x) for x in stopping_power_MeV_list]) #från MeV g/cm^2 till eV kg/m^2

    CSDA_g_per_cm2_list=stopping_power_data[:, 2]
    CSDA_list=np.array([(lambda x: x * 10**1)(x) for x in CSDA_g_per_cm2_list]) #från g/cm^2 till kg/m^2

    # Tar index för närmaste energi på alfapartikeln.
    diff = np.abs(energi_list - energi)
    closest_indices = np.argsort(diff)[:2]

    # Närmsta värdena:
    energi_close = energi_list[closest_indices]
    stopping_power_close = STP_list[closest_indices]
    CSDA_close = CSDA_list[closest_indices]

    # Linjärinterpolera och få fram stopping power och slumpmässig steglängd.
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        stopping_power = stopping_power_close[0]
        CSDA = CSDA_close[0]

    else:
        stopping_power = (stopping_power_close[0] + (energi - energi_close[0]) * (
                stopping_power_close[1] - stopping_power_close[0]) / (
                                  energi_close[1] - energi_close[0]))

        CSDA = (CSDA_close[0] + (energi - energi_close[0]) * (
                CSDA_close[1] - CSDA_close[0]) / (
                        energi_close[1] - energi_close[0]))

    steglängd = CSDA / rho_medium * 10 ** (-2)  # Omvandla cm till m.
    stopping_power = stopping_power * rho_medium * 10 ** 2  # Omvandla till eV / m.

    return stopping_power, steglängd
