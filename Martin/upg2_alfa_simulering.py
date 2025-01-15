from imports import *
from upg1_sampla_energi_start import energi_start
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd
from upg2_riktning import riktning_uniform, riktning_skal
from upg2_position_start import position_start_innanför, position_start_skal
from upg2_laddad_partikel_väg import laddad_partikel_väg


def run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_alpha, radie_sfär,
                 max_antal_steg):
    """
    Monte-Carlo simulering för alfapartiklarna.
    :param iterationer: Antal sönderfall som ska simuleras.
    :param rho_medium: Densiteten av mediumet (vatten).
    :param radie_partikel: Partikelns radie (alfa).
    :param stopping_power_data: Stopping power data.
    :param position_start_alpha: Fördelningsfunktion: uniform eller ytfördelning.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Summeringen av energideponeringen innanför sfären.
    """

    energideponering_summa = 0
    #   ----------------------------------------------------------------------
    #   Vilken fördelningsfunktion som ska användas bestämmer hur
    #   sampling av riktning och position sker.
    #   ----------------------------------------------------------------------
    if position_start_alpha == position_start_skal:
        #   ----------------------------------------------------------------------
        #   Ytfördelning på en sfär.
        #   ----------------------------------------------------------------------
        iterationer = 0.5 * iterationer
        for i in range(int(iterationer)):

            # Sampla startenergin.
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6)} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_skal()
            position_start = position_start_skal(radie_sfär, radie_partikel)

            # Sampla steglängd för partikeln.
            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering = laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, stopping_power_data, max_antal_steg)

            # Summera alla dosbidrag.
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    else:
        #   ----------------------------------------------------------------------
        #   Uniform fördelning i en sfär.
        #   ----------------------------------------------------------------------
        for i in range(iterationer):

            # Sampla startenergin.
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6)} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_uniform()
            position_start = position_start_innanför(radie_sfär)

            # Sampla steglängd för partikeln.
            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering = laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, stopping_power_data, max_antal_steg)

            # Summera alla dosbidrag.
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    # print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
    return energideponering_summa


def energideponering_eV_till_Gy(energideponering_eV, rho_medium, radie_sfär):
    """
    Beräkna absorberad i termer av Gy.
    :param energideponering_eV: Absorberad dos i eV.
    :param rho_medium: Mediumets densitet.
    :param radie_sfär: Sfärens radie.
    """

    # Beräkna volym (m^3) och massa (kg).
    V = 4 / 3 * np.pi * radie_sfär ** 3
    massa = V * rho_medium

    # Beräkna absorberad dos i Joule.
    energideponering_J = energideponering_eV * 1.602 * 10 ** (-19)

    # Beräkna absorberad dos i Gy.
    energideponering_Gy = energideponering_J / massa

    return energideponering_Gy


if __name__ == "__main__":
    #   ----------------------------------------------------------------------
    #   Kör simuleringen med ingångsvärden.
    #   ----------------------------------------------------------------------
    iterationer = 10 ** 2
    dummy_iterationer = 10 ** 1
    max_antal_steg = 10 ** 3

    stopping_power_data = np.loadtxt(stopping_power_alfa_file)

    rho_medium = rho_vatten
    radie_partikel = radie_alpha

    radie_sfär_skal = 300 * 10 ** (-6)
    radie_sfär_innanför = 1 * 10 ** (-3)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_skal,
                     radie_sfär_skal, max_antal_steg)

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                             position_start_skal, radie_sfär_skal, max_antal_steg)
    energideponering_skal_Gy = energideponering_eV_till_Gy(energideponering_tot_skal, rho_medium, radie_sfär_skal)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_innanför,
                     radie_sfär_innanför, max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_alpha(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                 position_start_innanför, radie_sfär_innanför, max_antal_steg)
    energideponering_innanför_Gy = energideponering_eV_till_Gy(energideponering_tot_innanför, rho_medium,
                                                               radie_sfär_innanför)
    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nRESULTAT\n----------------------------------------------------------------------\n')

    #   ----------------------------------------------------------------------
    #   Beräkna resultat och jämför med valideringsdata.
    #   ----------------------------------------------------------------------
    print(
        f'\nSkal (300 mikrometer): Energideponering:\n{energideponering_skal_Gy * 10 ** 6 / iterationer} E-06 Gy / sönderfall')
    print(f'faktor {(energideponering_skal_Gy * 10 ** 6 / iterationer) / 1.66:.3f} av facit')

    print(
        f'\nInnanför (1 mm): Energideponering:\n{energideponering_innanför_Gy * 10 ** 8 / iterationer} E-08 Gy / sönderfall')
    print(f'faktor {(energideponering_innanför_Gy * 10 ** 8 / iterationer) / 9.18:.3f} av facit')
