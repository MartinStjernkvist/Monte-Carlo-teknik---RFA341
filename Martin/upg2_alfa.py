from imports import *
from upg1_sampla_energi_start import energi_start
from MC_Linnea.Alfa_stp_och_RCSDA import Stopping_power_och_steglängd


@jit(nopython=True)
def riktning_alpha():
    """
    Uniform riktning.
    """
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi


@jit(nopython=True)
def riktning_alpha_skal():
    """
    Ifall skalförfördelnning:
    pi / 2 < phi < 3 * pi / 2 för att effektivisera koden.
    Då kan antalet iterationer halveras, eftersom man vet att
    hälften av partiklarna ändå skulle lämna sfären.
    """
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = pi / 4 * (np.random.rand() + 2)
    return theta, phi


"""
def steglängd_alpha(energi, df_stopping_power):
    # Avstånd till braggtopp
    # print('WIP')
    medelvägslängd = 10**(-6)
    return medelvägslängd
"""


@jit(nopython=True)
def energiförlust_alpha(energi, steg):
    """
    Funktion som beräknar energiförlusten för alfapartiklarna.
    """
    # Implementera stopping power
    STP, _ = Stopping_power_och_steglängd(energi)
    energiförlust = STP * steg  # i MeV
    energi -= energiförlust

    if energi <= 0:
        energi = 0

    return energi


@jit(nopython=True)
def position_start_alpha_innanför(radie_sfär, phi, theta):
    r = radie_sfär * np.random.rand()

    x = r * np.sin(theta) * np.cos(phi)
    # Utnyttja sfärisk symmetri.
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
def position_start_alpha_skal(radie_sfär, phi, theta):
    r = radie_sfär - 0.5 * radie_alpha  # För att inte endast theta = pi ska ge utslag

    x = r * np.sin(theta) * np.cos(phi)
    # Utnyttja sfärisk symmetri.
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
def laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie, max_antal_steg=100):
    """
    Funktion som följer alfapartikeln allteftersom den växelverkar i ett medium.
    :param radie: Radien av sfären för fördelningen.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Energideponeringen innanför sfären.
    """

    position_vektor = start_position
    energi = start_energi

    steg_storlek = steglängd / max_antal_steg

    riktning = np.array(
        [np.sin(theta) * np.cos(phi)
            , np.sin(theta) * np.cos(phi)
            , np.cos(theta)])

    riktning /= np.linalg.norm(riktning)
    steg_vektor = riktning * steg_storlek

    # Under tiden som partikeln fortfarande inte tagit hela sitt steg.
    for i in range(max_antal_steg):

        position_vektor += steg_vektor
        energi = energi - energiförlust_alpha(energi, steg_storlek)

        if np.dot(position_vektor, position_vektor) <= radie:
            print(f'Energideponering i position ', position_vektor)
        else:
            break

    energideponering = start_energi - energi

    return energideponering


def run_MC_alpha(iterationer, df_stopping_power, position_start_alpha, radie, max_antal_steg):
    """
    Monte-Carlo simulering för alfapartiklarna.
    :param iterationer: Antal sönderfall som ska simuleras.
    :param df_stopping_power: Stopping power data.
    :param position_start_alpha: Uniform fördelning i sfären, eller ytfördelning.
    :param radie: Radien av sfären för fördelningen.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Summeringen av energideponeringen innanför sfären.
    """

    energideponering_summa = 0
    start_energi = energi_start(At211_energi, At211_sannolikhet)

    if position_start_alpha == position_start_alpha_skal:
        iterationer = 0.5 * iterationer
        for i in range(iterationer):
            theta, phi = riktning_alpha_skal()
            start_position = position_start_alpha(radie, phi, theta)
            _, steglängd = Stopping_power_och_steglängd(
                start_energi)  # steglängd_alpha(start_position, df_stopping_power)
            energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                   max_antal_steg)

    else:
        for i in range(iterationer):
            theta, phi = riktning_alpha()
            start_position = position_start_alpha(radie, phi, theta)
            _, steglängd = Stopping_power_och_steglängd(
                start_energi)  # steglängd_alpha(start_position, df_stopping_power)
            energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                   max_antal_steg)

        energideponering_summa += energideponering

    print(f'\nEnergideponering per partikel: {energideponering_summa / iterationer:.2f} eV / partikel')
    return energideponering_summa


if __name__ == "__main__":
    iterationer = 10 ** 2
    dummy_iterationer = 10 ** 2
    max_antal_steg = 10 ** 3

    df_stopping_power = pd.read_excel(attenueringsdata_file)

    radie_sfär = 300 * 10 ** (-6)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, df_stopping_power, position_start_alpha_skal, radie_sfär,
                     max_antal_steg)

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_alpha(iterationer, df_stopping_power, position_start_alpha_skal, radie_sfär,
                                             max_antal_steg)

    end_time(start)

    radie_sfär = 1 * 10 ** (-3)

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, df_stopping_power, position_start_alpha_innanför, radie_sfär,
                     max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_alpha(iterationer, df_stopping_power, position_start_alpha_innanför,
                                                 radie_sfär,
                                                 max_antal_steg)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nRESULTAT\n----------------------------------------------------------------------\n')

    print(f'\nSkal: Energideponering per partikel: {energideponering_tot_skal / iterationer:.2f} eV / partikel')
    print(f'Innanför: Energideponering per partikel: {energideponering_tot_innanför / iterationer:.2f} eV / partikel')
