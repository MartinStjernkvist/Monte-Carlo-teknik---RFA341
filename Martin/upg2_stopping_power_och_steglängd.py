from imports import *


def stopping_power_och_steglängd(energi, rho, stopping_power_data):
    """
    Funktion som ger stopping power och
    :param energi: Partikelns energi i eV.
    :return: Stopping power (eV/m) och steglängden (m).
    """

    data_energi = stopping_power_data[:, 0]
    data_stoppingpower = stopping_power_data[:, 1]
    data_range = (stopping_power_data[:, 2])

    # Tar index för närmaste energi på alfapartikeln.
    diff = np.abs(data_energi - energi)
    closest_indices = np.argsort(diff)[:2]

    # Närmsta värdena:
    energi_close = data_energi[closest_indices]
    stopping_power_close = data_stoppingpower[closest_indices]
    range_close = data_range[closest_indices]

    # Linjärinterpolera och få fram stopping power och slumpmässig steglängd.
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        stopping_power = stopping_power_close[0] * rho
        steglängd = range_close[0] / rho

    else:
        stopping_power = (stopping_power_close[0] + (energi - energi_close[0]) * (
                stopping_power_close[1] - stopping_power_close[0]) / (
                                  energi_close[1] - energi_close[0])) * rho

        steglängd = (range_close[0] + (energi - energi_close[0]) * (
                range_close[1] - range_close[0]) / (
                             energi_close[1] - energi_close[0])) / rho

    steglängd = steglängd * 10 ** (-2)  # Omvandla till m.
    stopping_power = stopping_power * 10 ** 2  # Omvandla till eV / m.

    return stopping_power, steglängd
