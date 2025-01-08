from imports import *
from upg1_sampla_energi_start import energi_start
from MC_Linnea.Alfa_stp_och_RCSDA import Stopping_power_och_steglängd

@jit(nopython=True)
def riktning_alpha():
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi


def steglängd_alpha(energi, df_stopping_power):
    # Avstånd till braggtopp
    # print('WIP')
    medelvägslängd = 10**(-6)
    return medelvägslängd


@jit(nopython=True)
def förflyttning(position_vektor, steg_vektor):
    position_vektor += steg_vektor
    return position_vektor


@jit(nopython=True)
def energiförlust_alpha(energi, steg):
    # Implementera stopping power
    # print('WIP')

    energiförlust = energi * 0.1
    energi -= energiförlust

    if energi <= 0:
        energi = 0

    return energi


@jit(nopython=True)
def position_start_alpha_innanför(radie_sfär, phi, theta):
    r = radie_sfär * np.random.rand()

    x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
def position_start_alpha_skal(radie_sfär, phi, theta):
    radie_alpha=1.2*10**(-15)*4**(1/3) #radie i meter enligt Physics Handbook
    r = radie_sfär - 0.5 * radie_alpha  # För att inte endast theta = pi ska ge utslag

    x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
# @jit(nopython=True)
def laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie, max_antal_steg=100):
    position_vektor = start_position
    energi = start_energi

    # trajectory = [tuple(position_vektor)]

    steg_storlek = steglängd / max_antal_steg

    riktning = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.cos(phi), np.cos(theta)])
    riktning /= np.linalg.norm(riktning)
    steg_vektor = riktning * steg_storlek


    # Under tiden som partikeln fortfarnade inte tagit hela sitt steg.
    for i in range(max_antal_steg):

        # print('steg_vektor', steg_vektor)
        position_vektor += steg_vektor
        energi = energi - energiförlust_alpha(energi, steg_storlek)

        if np.dot(position_vektor, position_vektor) <= radie:
            innanför = True
            # trajectory.append(tuple(position_vektor))
            print(f'Energideponering i position ', position_vektor)
        else:
            break
            # print('Partikel utanför sfär!')

    energideponering = start_energi - energi

    return energideponering  # , trajectory


def run_MC_alpha(iterationer, df_stopping_power, position_start_alpha, radie, max_antal_steg):
    energideponering_summa = 0
    utanför = 0
    start_energi==energi_start(At211_energi,At211_sannolikhet)

    if position_start_alpha == position_start_alpha_skal:

        for i in range(iterationer):
            theta, phi = riktning_alpha()

            if not pi / 2 < phi < 3 * pi / 2:
                # print('Utanför')
                utanför += 1
                energideponering = 0
            else:
                start_position = position_start_alpha(radie, phi, theta)
                steglängd = steglängd_alpha(start_position, df_stopping_power)
                energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                       max_antal_steg)

    else:
        for i in range(iterationer):
            theta, phi = riktning_alpha()
            start_position = position_start_alpha(radie, phi, theta)
            steglängd = steglängd_alpha(start_position, df_stopping_power)
            energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                   max_antal_steg)

        energideponering_summa += energideponering

    print('antal utanför: ', utanför)
    print('total energideponering: ', energideponering_summa)
    print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
    return energideponering_summa

#
# def run_MC_alpha_innanför(iterationer, df_stopping_power, start_energi, radie, max_antal_steg):
#     energideponering_summa = 0
#     utanför = 0
#
#     for i in range(iterationer):
#         theta, phi = riktning_alpha()
#
#         if not pi / 2 < phi < 3 * pi / 2:
#             # print('Utanför')
#             utanför += 1
#             energideponering = 0
#         else:
#             start_position = position_start_alpha_innanför(radie, phi, theta)
#             steglängd = steglängd_alpha(start_position, df_stopping_power)
#             energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
#                                                    max_antal_steg)
#
#         energideponering_summa += energideponering
#
#     print('antal utanför: ', utanför)
#     print('total energideponering: ', energideponering_summa)
#     print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
#     return energideponering_summa


if __name__ == "__main__":
    iterationer = 10 ** 2
    dummy_iterationer = 10**2
    max_antal_steg = 10**3

    df_stopping_power = pd.read_excel(attenueringsdata_file)

    radie_sfär = 300 * 10 ** (-6)
    start_energi = 10**4

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, df_stopping_power, position_start_alpha_skal, start_energi, radie_sfär,
                                        max_antal_steg)

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_alpha(iterationer, df_stopping_power, position_start_alpha_skal, start_energi, radie_sfär, max_antal_steg)

    end_time(start)

    radie_sfär = 1 * 10 ** (-3)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, df_stopping_power, position_start_alpha_innanför, start_energi, radie_sfär,
                                                 max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_alpha(iterationer, df_stopping_power, position_start_alpha_innanför, start_energi, radie_sfär,
                                                 max_antal_steg)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nRESULTAT\n----------------------------------------------------------------------\n')

    print(f'\nSkal: Energideponering per partikel: {energideponering_tot_skal / iterationer:.2f} eV / partikel')
    print(f'Innanför: Energideponering per partikel: {energideponering_tot_innanför / iterationer:.2f} eV / partikel')

