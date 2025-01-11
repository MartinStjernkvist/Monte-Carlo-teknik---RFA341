from imports import *
from upg1_sampla_energi_start import energi_start
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd
from upg2_riktning import riktning_uniform, riktning_skal
from upg2_position_start import position_start_innanför, position_start_skal
from upg2_energi_efter_förlust import energi_efter_energiförlust
from upg2_laddad_partikel_väg import laddad_partikel_väg


def run_MC_alpha(iterationer, rho_medium, radie_partikel, df_stopping_power, position_start_alpha, radie_sfär,
                 max_antal_steg):
    """
    Monte-Carlo simulering för alfapartiklarna.
    :param iterationer: Antal sönderfall som ska simuleras.
    :param df_stopping_power: Stopping power data.
    :param position_start_alpha: Uniform fördelning i sfären, eller ytfördelning.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Summeringen av energideponeringen innanför sfären.
    """

"""
    energideponering_summa = 0
    x_lista=[]
    y_lista=[]
    z_lista=[]
    dos_lista=[]

    if position_start_alpha == position_start_skal:
        iterationer = 0.5 * iterationer
        for i in range(int(iterationer)):
            energi = energi_start(At211_energi, At211_sannolikhet)

            theta, phi = riktning_skal()
            position_start = position_start_alpha(radie_sfär, radie_partikel)

            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, df_stopping_power)

            energideponering= laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, df_stopping_power, max_antal_steg)
            
            x_lista+=x
            y_lista+=y
            z_lista+=z
            dos_lista+=dos
            
            energideponering_summa += energideponering
    else:
        for i in range(iterationer):
            energi = energi_start(At211_energi, At211_sannolikhet)

            theta, phi = riktning_uniform()
            position_start = position_start_alpha(radie_sfär)

            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, df_stopping_power)

            energideponering = laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, df_stopping_power, max_antal_steg)
            #, x,y,z,dos
            
            x_lista+=x
            y_lista+=y
            z_lista+=z
            dos_lista+=dos
           
            energideponering_summa += energideponering

    # print('antal utanför: ', utanför)

    #Plotta en figur som visar energideponeringen i hela tumören
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x_lista,y_lista,z_lista,c=dos_lista,cmap='plasma')
    #plt.show()
    #plt.colorbar()
    print('total energideponering: ', energideponering_summa)
    print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
    return energideponering_summa


def energideponering_eV_till_Gy(energideponering_eV, rho_medium, radie_sfär):
    V = 4 / 3 * np.pi * radie_sfär ** 3
    rho_kg_m3 = rho_medium * 1000

    massa = V * rho_kg_m3
    energideponering_J = energideponering_eV * 1.602 * 10 ** (-19)
    energideponering_Gy = energideponering_J / massa  # J / kg

    return energideponering_Gy


if __name__ == "__main__":
    iterationer = 10 ** 4
    dummy_iterationer = 10 ** 2
    max_antal_steg = 10 ** 4

    stopping_power_data = np.loadtxt('MC_Linnea/Stoppingpower_data_alfa')

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
    energideponering_innanför_Gy = energideponering_eV_till_Gy(energideponering_tot_innanför, rho_medium, radie_sfär_innanför)
    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nRESULTAT\n----------------------------------------------------------------------\n')

    print(f'\nSkal (300 mikrometer): Energideponering: {energideponering_skal_Gy * 10 ** 9/ iterationer:.2f} E-09  Gy / sönderfall')
    print(f'Innanför (1 mm): Energideponering: {energideponering_tot_innanför* 10 ** 9 / iterationer:.2f} E-09 Gy / sönderfall')

    """
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
    :param stopping_power_data: Stopping power data.
    :param position_start_alpha: Uniform fördelning i sfären, eller ytfördelning.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Summeringen av energideponeringen innanför sfären.
    """

    energideponering_summa = 0

    if position_start_alpha == position_start_skal:
        iterationer = 0.5 * iterationer
        for i in range(int(iterationer)):
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6)} MeV')

            theta, phi = riktning_skal()
            position_start = position_start_skal(radie_sfär, radie_partikel)

            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            print(f'steglängd: {steglängd * 10 ** (-6):.2f} mikrometer')

            energideponering = laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, stopping_power_data, max_antal_steg)

            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    else:
        for i in range(iterationer):
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6)} MeV')

            theta, phi = riktning_uniform()
            position_start = position_start_innanför(radie_sfär)

            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            energideponering = laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, stopping_power_data, max_antal_steg)

            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    # print('antal utanför: ', utanför)
    # print('total energideponering: ', energideponering_summa)
    # print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
    return energideponering_summa


def energideponering_eV_till_Gy(energideponering_eV, rho_medium, radie_sfär):
    V = 4 / 3 * np.pi * radie_sfär ** 3
    massa = V * rho_medium  # kg
    energideponering_J = energideponering_eV * 1.602 * 10 ** (-19)  # J

    energideponering_Gy = energideponering_J / massa  # J / kg

    return energideponering_Gy


if __name__ == "__main__":
    iterationer = 10 ** 2
    dummy_iterationer = 10 ** 1
    max_antal_steg = 10 ** 3

    stopping_power_data = np.loadtxt('MC_Linnea/Stoppingpower_data_alfa')

    rho_medium = rho_vatten
    radie_partikel = radie_alpha

    radie_sfär_skal = 300 * 10 ** (-6) #m
    radie_sfär_innanför = 1 * 10 ** (-3)#m

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

    print(
        f'\nSkal (300 mikrometer): Energideponering:\n{energideponering_skal_Gy * 10 ** 6 / iterationer} E-06 Gy / sönderfall')
    print(f'faktor {(energideponering_skal_Gy * 10 ** 6 / iterationer) / 1.66} av facit')

    print(
        f'\nInnanför (1 mm): Energideponering:\n{energideponering_innanför_Gy * 10 ** 8 / iterationer} E-08 Gy / sönderfall')
    print(f'faktor {(energideponering_innanför_Gy * 10 ** 8 / iterationer) / 9.18} av facit')
