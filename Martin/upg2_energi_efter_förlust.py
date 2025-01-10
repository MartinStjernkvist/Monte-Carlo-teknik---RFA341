from imports import *
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd

def energi_efter_energiförlust(energi, steglängd, rho, stopping_power_data):
    """
    Funktion som beräknar energiförlusten efter varje steg som partikeln tar.
    :param energi: Partikelns energi innan steget.
    :param steglängd: Stegstorleken (cm).
    :param rho: Densiteten av mediumet.
    :param stopping_power_data: Stopping power data.
    :return: Energiförlusten under steget.
    """
    stopping_power, _ = stopping_power_och_steglängd(energi, rho, stopping_power_data)
    energiförlust = stopping_power * steglängd  # i eV
    ny_energi = energi - energiförlust

    if ny_energi <= 0:
        ny_energi = 0

    return ny_energi