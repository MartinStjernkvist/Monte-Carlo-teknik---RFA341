from imports import *
from upg2_elektron_energi import elektron_startenergi, elektron_energi_start
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd
from upg2_riktning import riktning_uniform, riktning_skal
from upg2_position_start import position_start_innanför, position_start_skal
from upg2_laddad_partikel_väg_elektron import laddad_partikel_väg_elektron
from upg2_alfa import energideponering_eV_till_Gy


def run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                    position_start_alpha, radie_sfär,
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
    x_list, y_list, z_list, dos_list = [], [], [], []

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
            # energi = elektron_startenergi() * 10 ** 6  # i eV
            energi = elektron_energi_start() * 10 ** 6  # i eV

            print(f'energi: {energi * 10 ** (-6)} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_skal()
            position_start = position_start_skal(radie_sfär, radie_partikel)

            # Sampla steglängd för partikeln.
            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            steglängd = steglängd * np.random.rand()
            # print(f'steglängd: {steglängd * 10 ** (-6):.2f} mikrometer')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            # kommer gå i Tau i riktningen men sen ändra på den i Steglängd-Tau med theta riktning
            energideponering, x, y, z, dos = laddad_partikel_väg_elektron(energi, position_start, phi, theta, steglängd,
                                                                          radie_sfär,
                                                                          rho_medium, stopping_power_data,
                                                                          scatter_power_data,
                                                                          max_antal_steg)
            # ändra i laddad_partikel_väg så den ändrar koordinatsystem med polarvinkeln!!!!!!!!!

            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    else:
        #   ----------------------------------------------------------------------
        #   Uniform fördelning i en sfär.
        #   ----------------------------------------------------------------------
        for i in range(iterationer):
            # Sampla startenergin.
            # energi = elektron_startenergi() * 10 ** 6  # i eV
            # energi = elektron_startenergi() * 10 ** 6  # i eV
            energi = elektron_energi_start() * 10 ** 6  # i eV
            print(f'energi: {energi * 10 ** (-6)} MeV')

            # Sampla riktning och startposition.
            theta, phi = riktning_uniform()
            position_start = position_start_innanför(radie_sfär)

            # Sampla steglängd för partikeln.
            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            steglängd = steglängd * np.random.rand()
            # print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            # Beräkna den totala energideponeringen för en partikel som växelverkar i sfären.
            energideponering, x, y, z, dos = laddad_partikel_väg_elektron(energi, position_start, phi, theta, steglängd,
                                                                          radie_sfär,
                                                                          rho_medium, stopping_power_data,
                                                                          scatter_power_data,
                                                                          max_antal_steg)

            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    # print('antal utanför: ', utanför)
    # print('total energideponering: ', energideponering_summa)
    # print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')

    """
    # Plotta en figur som visar energideponeringen i hela tumören
    # Hitta ett sätt att färga punkterna för värde på energideponeringen
    fig = plt.figure(1)
    fig2 = plt.figure(2)
    ax = fig.add_subplot(projection='3d')
    ax2 = fig2.add_subplot(projection='3d')
    ax.scatter(x_list, y_list, z_list, c=dos_list, cmap='plasma')
    ax.set_xlabel('x-axel i meter')
    ax.set_ylabel('y-axel i meter')

    # Testar att sätta en sfär för tumören
    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = radie_sfär * np.cos(u) * np.sin(v)
    y = radie_sfär * np.sin(u) * np.sin(v)
    z = radie_sfär * np.cos(v)
    ax.plot_wireframe(x, y, z, color="k", alpha=0.3)

    # ax2.scatter(x_list,y_list,dos_list,c=dos_list,cmap='plasma')
    ax2.hist(dos_list)
    ax2.set_xlabel('x-axel i energideponering i eV')
    ax2.set_ylabel('y-axel i antal som har samma energi')
    plt.show()
    """
    return energideponering_summa


if __name__ == "__main__":
    iterationer = 10 ** 4
    dummy_iterationer = 10 ** 2
    max_antal_steg = 10 ** 3

    stopping_power_data = np.loadtxt(elektron_stopping_power_data)
    scatter_power_data = np.loadtxt(elektron_scatter_power_data)

    rho_medium = rho_vatten
    radie_partikel = r_e_m

    radie_sfär_skal = 300 * 10 ** (-6)  # m
    radie_sfär_innanför = 1 * 10 ** (-3)  # m

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_elektron(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                        position_start_skal,
                        radie_sfär_skal, max_antal_steg)

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                scatter_power_data,
                                                position_start_skal, radie_sfär_skal, max_antal_steg)
    energideponering_skal_Gy = energideponering_eV_till_Gy(energideponering_tot_skal, rho_medium, radie_sfär_skal)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_elektron(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, scatter_power_data,
                        position_start_innanför,
                        radie_sfär_innanför, max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_elektron(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                                    scatter_power_data,
                                                    position_start_innanför, radie_sfär_innanför, max_antal_steg)
    energideponering_innanför_Gy = energideponering_eV_till_Gy(energideponering_tot_innanför, rho_medium,
                                                               radie_sfär_innanför)
    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nRESULTAT\n----------------------------------------------------------------------\n')

    print(
        f'\nSkal (300 mikrometer): Energideponering:\n{energideponering_skal_Gy * 10 ** 8 / iterationer} E-08 Gy / sönderfall')
    print(f'faktor {(energideponering_skal_Gy * 10 ** 8 / iterationer) / 4.07} av facit')

    print(
        f'\nInnanför (1 mm): Energideponering:\n{energideponering_innanför_Gy * 10 ** 9 / iterationer} E-09 Gy / sönderfall')
    print(f'faktor {(energideponering_innanför_Gy * 10 ** 9 / iterationer) / 5.22} av facit')
