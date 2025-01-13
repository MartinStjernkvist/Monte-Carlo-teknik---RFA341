from imports import *
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd


def energi_efter_energiförlust(energi, steglängd, rho_medium, stopping_power_data):
    """
    Funktion som beräknar energiförlusten efter varje steg som partikeln tar.
    :param energi: Partikelns energi innan steget.
    :param steglängd: Stegstorleken (m).
    :param rho_medium: Densiteten av mediumet (kg/m^3).
    :param stopping_power_data: Stopping power data.
    :return: Energiförlusten under steget.
    """

    # Beräkna stopping power, utifrån tabellerad data.
    stopping_power, _ = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
    # print(f'stopping power: {stopping_power}')

    # Beräkna energiförlusten med hjälp av stopping power.
    energiförlust = stopping_power * steglängd
    # print(f'energiförlust: {energiförlust * 10**(-6):.4f} MeV')

    # Beräkna partikelns nya energi, efter steget har tagits.
    ny_energi = energi - energiförlust

    # Om partikelns energi < 0, ansätt att partikelns energi = 0.
    if ny_energi <= 0:
        ny_energi = 0

    return ny_energi
