from imports import *
from Elektron_stp_och_steglängd import Stopping_power_och_steglängd_elektron

def energi_efter_energiförlust(energi, steglängd, rho_medium, stopping_power_data):
    """
    Funktion som beräknar energiförlusten efter varje steg som partikeln tar.
    :param energi: Partikelns energi innan steget.
    :param steglängd: Stegstorleken (cm).
    :param rho_medium: Densiteten av mediumet.
    :param stopping_power_data: Stopping power data.
    :return: Energiförlusten under steget.
    """
    stopping_power, _ ,_= Stopping_power_och_steglängd_elektron(energi, rho_medium, stopping_power_data)
    energiförlust = stopping_power * steglängd  # i eV
    ny_energi = energi - energiförlust

    if ny_energi <= 0:
        ny_energi = 0

    return ny_energi