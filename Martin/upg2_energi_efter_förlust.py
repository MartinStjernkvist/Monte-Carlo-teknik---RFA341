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

    stopping_power, _ = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
    print(f'stopping power: {stopping_power}')

    energiförlust = stopping_power * steglängd
    print(f'energiförlust: {energiförlust * 10**(-6):.2f} MeV')
    ny_energi = energi - energiförlust

    if ny_energi <= 0:
        ny_energi = 0

    return ny_energi