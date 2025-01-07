from imports import *


@jit(nopython=True)
def riktning_alpha():
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi


def steglängd_alpha(energi, df_stopping_power):
    # Avstånd till braggtopp
    # print('WIP')
    medelvägslängd = 100
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
    r = radie_sfär - 0.5 * radie_alpha  # För att inte endast theta = pi ska ge utslag

    x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
def position_start_alpha_skal(radie_sfär, phi, theta):
    r = radie_sfär * np.random.rand()

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

    innanför = True
    summa_steglängd = 0

    # Under tiden som partikeln fortfarnade inte tagit hela sitt steg.
    for i in range(max_antal_steg):

        # Medan partikeln fortfarande är innanför tumörcellen.
        while innanför == True:

            # Det sista värdet kommer sparas, resterande är irrelevanta.

            # print('steg_vektor', steg_vektor)
            position_vektor += steg_vektor
            energi = energi - energiförlust_alpha(energi, steg_storlek)

            if np.dot(position_vektor, position_vektor) <= radie:
                innanför = True
                # trajectory.append(tuple(position_vektor))
                x, y, z = position_vektor[0], position_vektor[1], position_vektor[2]
                print(f'Energideponering i position [{x, y, z}].')
            else:
                innanför = False
                # print('Partikel utanför sfär!')

    energideponering = start_energi - energi

    return energideponering  # , trajectory


def run_MC_alpha_skal(iterationer, df_stopping_power, start_energi, radie, max_antal_steg):
    energideponering_summa = 0
    utanför = 0

    for i in range(iterationer):
        theta, phi = riktning_alpha()

        if not pi / 2 < phi < 3 * pi / 2:
            # print('Utanför')
            utanför += 1
            energideponering = 0
        else:
            start_position = position_start_alpha_skal(radie, phi, theta)
            steglängd = steglängd_alpha(start_position, df_stopping_power)
            energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                   max_antal_steg)

        energideponering_summa += energideponering

    print('antal utanför: ', utanför)
    print('total energideponering: ', energideponering_summa)
    print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
    return energideponering_summa


def run_MC_alpha_innanför(iterationer, df_stopping_power, start_energi, radie, max_antal_steg):
    energideponering_summa = 0
    utanför = 0

    for i in range(iterationer):
        theta, phi = riktning_alpha()

        if not pi / 2 < phi < 3 * pi / 2:
            # print('Utanför')
            utanför += 1
            energideponering = 0
        else:
            start_position = position_start_alpha_innanför(radie, phi, theta)
            steglängd = steglängd_alpha(start_position, df_stopping_power)
            energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                   max_antal_steg)

        energideponering_summa += energideponering

    print('antal utanför: ', utanför)
    print('total energideponering: ', energideponering_summa)
    print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
    return energideponering_summa


if __name__ == "__main__":
    iterationer = 10 ** 5
    dummy_iterationer = 10**3
    max_antal_steg = 10**5

    df_stopping_power = pd.read_excel(attenueringsdata_file)

    radie_sfär = 300 * 10 ** (-6)  # byt
    start_energi = 1_000_000

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    energideponering_tot = run_MC_alpha_skal(dummy_iterationer, df_stopping_power, start_energi, radie_sfär,
                                             max_antal_steg)

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot = run_MC_alpha_skal(iterationer, df_stopping_power, start_energi, radie_sfär, max_antal_steg)

    end_time(start)

    radie_sfär = 10 ** (-3)  # byt
    start_energi = 1_000_000

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    energideponering_tot = run_MC_alpha_innanför(dummy_iterationer, df_stopping_power, start_energi, radie_sfär,
                                                 max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot = run_MC_alpha_innanför(iterationer, df_stopping_power, start_energi, radie_sfär,
                                                 max_antal_steg)

    end_time(start)
