from imports import *
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd


def steglängd_från_energi(energi_start, rho_medium, stopping_power_data, energiförlust_faktor):
    """
    Funktion som beräknar steglängd utifrån en energiförlust.
    """

    energi_1 = energi_start

    # Bestäm ny energi, varefter steglängden kan lösas ut m.h.a. stopping power.
    energi_2 = energi_1 * energiförlust_faktor  # antag energiförlust innan nästa kollission

    stopping_power_1, _ = stopping_power_och_steglängd(energi_1, rho_medium, stopping_power_data)
    stopping_power_2, _ = stopping_power_och_steglängd(energi_2, rho_medium, stopping_power_data)

    stopping_power_medel = (stopping_power_2 + stopping_power_1) / 2
    print(f'stopping power medel: {stopping_power_medel}')

    steglängd = - (energi_2 - energi_1) / stopping_power_medel
    print(f'steglängd: {steglängd * 10**6:.2f} mikrometer')

    return steglängd