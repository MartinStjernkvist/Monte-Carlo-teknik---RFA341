from imports import *
from upg2_elektron_energi import elektron_energi_start
from upg2_riktning import riktning_uniform, riktning_skal
from upg2_position_start import position_start_innanför, position_start_skal
from upg2_fixad_laddad_partikel_väg_elektroner import laddad_partikel_väg_elektron


def run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                    position_start, radie_sfär, energiförlust_faktor):
    """
    Monte-Carlo simulering för elektrontransport.
    :param iterationer: Antal sönderfall som ska simuleras.
    :param rho_medium: Mediumets densitet (kg / m^3).
    :param radie_partikel: Partikelns radie.
    :param stopping_power_data: Stopping power tabelldata.
    :param scatter_power_data: Scatter power tabelldata.
    :param position_start: Uniform fördelning i sfären, eller ytfördelning.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param energiförlust_faktor: Energiförlustfaktor efter varje steglängd.
    :return: Summeringen av energideponeringen innanför sfären.
    """

    # Initiera en energisumma och tomma listor för att spara datan.
    energideponering_summa = 0
    x_list, y_list, z_list, dos_list = [], [], [], []

    #   -----------------------------------
    #   Vilken fördelningsfunktion som ska användas bestämmer hur
    #   sampling av riktning och position sker.
    #   -----------------------------------
    if position_start == position_start_skal:
        #   -----------------------------------
        #   Ytfördelning på en sfär.
        #   -----------------------------------
        iterationer = 0.5 * iterationer
        for i in range(int(iterationer)):
            # Sampla riktning och startposition.
            theta, phi = riktning_skal()
            position_start = position_start_skal(radie_sfär, radie_partikel)

            # Sampla startenergin.
            energi_start = elektron_energi_start() * 10 ** 6  # i eV

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering, x, y, z, dos = laddad_partikel_väg_elektron(energi_start, position_start, phi, theta,
                                                                          radie_sfär, rho_medium,
                                                                          stopping_power_data, scatter_power_data,
                                                                          energiförlust_faktor)

            # Spara mätpunkter för plottning.
            x_list += x
            y_list += y
            z_list += z
            dos_list += dos

            # Summera energideponeringsbidragen från respektive iteration.
            energideponering_summa += energideponering
            # print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    else:
        #   -----------------------------------
        #   Uniform fördelning i en sfär.
        #   -----------------------------------
        for i in range(iterationer):
            # Sampla startenergin.
            energi_start = elektron_energi_start() * 10 ** 6  # i eV

            # Sampla riktning och startposition.
            theta, phi = riktning_uniform()
            position_start = position_start_innanför(radie_sfär)

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering, x, y, z, dos = laddad_partikel_väg_elektron(energi_start, position_start, phi, theta,
                                                                          radie_sfär, rho_medium,
                                                                          stopping_power_data, scatter_power_data,
                                                                          energiförlust_faktor)
            # Spara mätpunkter för plottning.
            x_list += x
            y_list += y
            z_list += z
            dos_list += dos

            # Summera energideponeringsbidragen från respektive iteration.
            energideponering_summa += energideponering
            # print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    # print('total energideponering: ', energideponering_summa)

    #   -----------------------------------
    #   Visualisera resultat i figur.
    #   -----------------------------------
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x_list, y_list, z_list, c=dos_list, cmap='plasma', label='Partikel position')
    # Fixa colorbar för att se energideponeringen i figuren

    # fig.colorbar(ax=ax, label='Energideponering',)

    ax.set_xlabel('x-axel (m)')
    ax.set_ylabel('y-axel (m)')
    ax.set_zlabel('z-axel (m)')

    # Testar att sätta en sfär för tumören
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = radie_sfär * np.cos(u) * np.sin(v)
    y = radie_sfär * np.sin(u) * np.sin(v)
    z = radie_sfär * np.cos(v)
    ax.plot_wireframe(x, y, z, color="k", alpha=0.3, label='Tumören')
    ax.legend()

    # Visa figur
    plt.tight_layout()
    plt.show()

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
    #   -----------------------------------
    #   Köra simuleringskoden.
    #   -----------------------------------
    iterationer = 10 ** 3
    dummy_iterationer = 3 * 10**2

    energiförlust_faktor = 0.96  # Energi efter en kollision = 96% av föregående energin.

    stopping_power_data = np.loadtxt(elektron_stopping_power_data)
    scatter_power_data = np.loadtxt(elektron_scatter_power_data)

    rho_medium = rho_vatten
    radie_partikel = r_e_m

    radie_sfär_skal = 300 * 10 ** (-6)  # m
    radie_sfär_innanför = 1 * 10 ** (-3)  # m

    print(
        '\n-----------------------------------\nDUMMY\n-----------------------------------\n')

    _ = run_MC_elektron(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                        position_start_skal,
                        radie_sfär_skal, energiförlust_faktor)

    start = time.time()

    print(
        '\n-----------------------------------\nRIKTIG\n-----------------------------------\n')
    energideponering_tot_skal = run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                scatter_power_data,
                                                position_start_skal,
                                                radie_sfär_skal, energiförlust_faktor)

    energideponering_skal_Gy = energideponering_eV_till_Gy(energideponering_tot_skal, rho_medium, radie_sfär_skal)

    end_time(start)

    print(
        '\n-----------------------------------\nDUMMY\n-----------------------------------\n')

    _ = run_MC_elektron(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                        position_start_innanför,
                        radie_sfär_innanför, energiförlust_faktor)

    start = time.time()
    print(
        '\n-----------------------------------\nRIKTIG\n-----------------------------------\n')
    energideponering_tot_innanför = run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                    scatter_power_data,
                                                    position_start_innanför, radie_sfär_innanför, energiförlust_faktor)

    energideponering_innanför_Gy = energideponering_eV_till_Gy(energideponering_tot_innanför, rho_medium,
                                                               radie_sfär_innanför)
    end_time(start)

    print(
        '\n-----------------------------------\nRESULTAT\n-----------------------------------\n')

    print(
        f'\nSkal (300 mikrometer): Energideponering:\n{energideponering_skal_Gy * 10 ** 8 / iterationer} E-08 Gy / sönderfall')
    print(f'faktor {(energideponering_skal_Gy * 10 ** 8 / iterationer) / 4.07} av facit')

    print(
        f'\nInnanför (1 mm): Energideponering:\n{energideponering_innanför_Gy * 10 ** 9 / iterationer} E-09 Gy / sönderfall')
    print(f'faktor {(energideponering_innanför_Gy * 10 ** 9 / iterationer) / 5.22} av facit')
