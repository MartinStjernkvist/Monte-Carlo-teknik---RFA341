
from imports import *
from upg1_sampla_energi_start import energi_start
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd
from upg2_riktning import riktning_uniform, riktning_skal
from upg2_position_start import position_start_innanför, position_start_skal
from upg2_laddad_partikel_väg import laddad_partikel_väg


def run_MC_alpha_2(iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_alpha, radie_sfär,
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
    x_list,y_list,z_list,dos_list=[],[],[],[]
    if position_start_alpha == position_start_skal:
        iterationer = 0.5 * iterationer
        for i in range(int(iterationer)):
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6)} MeV')

            theta, phi = riktning_skal()
            position_start = position_start_skal(radie_sfär, radie_partikel)

            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            #print(f'steglängd: {steglängd * 10 ** (-6):.2f} mikrometer')

            energideponering , x,y,z, dos= laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, stopping_power_data, max_antal_steg)
            x_list+=x
            y_list+=y
            z_list+=z
            dos_list+=dos
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

    else:
        for i in range(iterationer):
            energi = energi_start(At211_energi, At211_sannolikhet)
            print(f'energi: {energi * 10 ** (-6)} MeV')

            theta, phi = riktning_uniform()
            position_start = position_start_innanför(radie_sfär)

            _, steglängd = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
            #print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

            energideponering , x,y,z, dos= laddad_partikel_väg(energi, position_start, phi, theta, steglängd, radie_sfär,
                                                   rho_medium, stopping_power_data, max_antal_steg)
            x_list+=x
            y_list+=y
            z_list+=z
            dos_list+=dos
            energideponering_summa += energideponering
            print(f'energideponering: {energideponering * 10 ** (-6)} MeV')

 

    # print('antal utanför: ', utanför)
    # print('total energideponering: ', energideponering_summa)
    # print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')

    #Plotta en figur som visar energideponeringen i hela tumören
    #Hitta ett sätt att färga punkterna för värde på energideponeringen
    fig = plt.figure(1)
    
    ax = fig.add_subplot(projection='3d')
    
    ax.scatter(x_list,y_list,z_list,c=dos_list,cmap='plasma',label='Partikel position')
    #Fixa colorbar för att se energideponeringen i figuren

    #fig.colorbar(ax=ax, label='Energideponering') 

    ax.set_xlabel('x-axel i meter')
    ax.set_ylabel('y-axel i meter')

    #Testar att sätta en sfär för tumören
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = radie_sfär*np.cos(u)*np.sin(v)
    y = radie_sfär*np.sin(u)*np.sin(v)
    z = radie_sfär*np.cos(v)
    ax.plot_wireframe(x, y, z, color="k",alpha=0.3,label='Tumören')
    ax.legend()
    

    #Visar dosfördelningen
    fig2=plt.figure(2)
    ax2=fig2.add_subplot()
    ax2.hist(dos_list)
    ax2.set_xlabel('Energideponering i eV')
    ax2.set_ylabel('Antal som har samma energi')

    #Visa figur
    plt.show()
    
    return energideponering_summa


def energideponering_eV_till_Gy(energideponering_eV, rho_medium, radie_sfär):
    V = 4 / 3 * np.pi * radie_sfär ** 3
    massa = V * rho_medium  # kg
    energideponering_J = energideponering_eV * 1.602 * 10 ** (-19)  # J

    energideponering_Gy = energideponering_J / massa  # J / kg

    return energideponering_Gy



if __name__ == "__main__":
    iterationer = 10 ** 3
    dummy_iterationer = 10 ** 2
    max_antal_steg = 10 ** 3

    stopping_power_data = np.loadtxt('MC_Linnea/Stoppingpower_data_alfa')

    rho_medium = rho_vatten
    radie_partikel = radie_alpha

    radie_sfär_skal = 300 * 10 ** (-6) #m
    radie_sfär_innanför = 1 * 10 ** (-3)#m

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha_2(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_skal,
                     radie_sfär_skal, max_antal_steg)

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_alpha_2(iterationer, rho_medium, radie_partikel, stopping_power_data,
                                             position_start_skal, radie_sfär_skal, max_antal_steg)
    energideponering_skal_Gy = energideponering_eV_till_Gy(energideponering_tot_skal, rho_medium, radie_sfär_skal)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha_2(dummy_iterationer, rho_medium, radie_partikel, stopping_power_data, position_start_innanför,
                     radie_sfär_innanför, max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_alpha_2(iterationer, rho_medium, radie_partikel, stopping_power_data,
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
