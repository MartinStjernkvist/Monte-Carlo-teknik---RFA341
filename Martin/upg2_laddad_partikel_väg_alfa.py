from imports import *
from upg2_energi_efter_förlust import energi_efter_energiförlust


def laddad_partikel_väg(energi_start, position_start, phi, theta, steglängd, radie_sfär, rho_medium,
                        stopping_power_data,
                        max_antal_steg):
    """
    Funktion som följer alfapartikeln allteftersom den växelverkar i ett medium.
    :param energi_start: Startenergi.
    :param position_start: Startposition.
    :param phi: Sfärisk vinkel.
    :param theta: Sfärisk vinkel.
    :param steglängd: Steglängd för partikeln i mediumet.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param rho_medium: Mediumets densitet.
    :param stopping_power_data: Tabellerad data för stopping power.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Energideponeringen innanför sfären.
    """

    # Initiera tomma listor för att spara datan.
    x, y, z, dos = [], [], [], []

    position_vektor = position_start
    energi = energi_start

    # Dela upp steglängden i mindre steg.
    steg_storlek = steglängd / max_antal_steg

    # Skapa en riktningsvektor.
    riktning = np.array(
        [np.sin(theta) * np.cos(phi)
            , np.sin(theta) * np.sin(phi)
            , np.cos(theta)])

    riktning /= np.linalg.norm(riktning)

    # Vektor för de små stegen som ska tas.
    steg_vektor = riktning * steg_storlek

    # Håll reda på antalet steg.
    steg_tagna = 0

    # Under tiden som partikeln fortfarande inte tagit hela sitt steg,
    # samt att partikeln inte förlorat all sin energi.
    while steg_tagna < max_antal_steg and energi > 0:

        steg_tagna += 1

        # Ta ett litet steg.
        position_vektor += steg_vektor

        # Håll reda på ifall partikeln befinner sig i sfären eller inte.
        # Ekvationen för en sfär: x^2 + y^2 + z^2 = r^2
        if np.sqrt(np.dot(position_vektor, position_vektor)) > radie_sfär:
            # print(f'Energideponering i position ', position_vektor)
            break

        else:
            # Beräkna energin efter steget.
            energi_ny = energi_efter_energiförlust(energi, steg_storlek, rho_medium, stopping_power_data)
            # print(f'ny energi: {energi * 10 ** (-6):.2f} MeV')

            # Spara mätpunkter.
            dos.append(energi - energi_ny)
            x.append(position_vektor[0])
            y.append(position_vektor[1])
            z.append(position_vektor[2])

            energi = energi_ny
            continue

    # Beräkna den totala energideponeringen (i sfären) från partikeln.
    energideponering = energi_start - energi
    # print(f'energideponering: {energideponering} eV')
    return energideponering, x, y, z, dos
