from imports import *
from upg2_energi_efter_förlust import energi_efter_energiförlust
from upg2_steg_transformation import ny_steg_transformera_koordinatsystem_3d
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd
from upg2_scattering_power import scattering_power_från_energi


def laddad_partikel_väg_elektron(energi_start, position_start, phi, theta, steglängd, radie_sfär, rho_medium,
                                 stopping_power_data, scatter_power_data,
                                 max_antal_steg):
    """
    Funktion som följer alfapartikeln allteftersom den växelverkar i ett medium.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Energideponeringen innanför sfären.
    """

    position_vektor = position_start
    energi = energi_start
    steg_tagna = 0
    x, y, z, dos = [], [], [], []

    # riktning = np.array(
    #     [np.sin(theta) * np.cos(phi)
    #         , np.sin(theta) * np.sin(phi)
    #         , np.cos(theta)])
    #
    # riktning /= np.linalg.norm(riktning)
    # position_vektor += riktning * steglängd


    # Under tiden som partikeln fortfarande inte tagit hela sitt steg.

    while steg_tagna < max_antal_steg and energi > 0:

        # Tar ett steg
        steg_tagna += 1

        # Nya värden på vinklarna
        phi_ny = np.random.random() * 2 * pi
        theta_ny = scattering_power_från_energi(energi, scatter_power_data, rho_medium)
        # stegstorlek totalt blir steglängd+cos(theta)*(s-steglängd) där s är CSDA
        _, s = stopping_power_och_steglängd(energi, rho_medium, stopping_power_data)
        tau = s * np.random.random()
        # tau = s * np.random.random() * 10**(-1)


        # Ändrar på positionsvektor efter att transformations matrisen
        dx, dy, dz = ny_steg_transformera_koordinatsystem_3d(steglängd, phi, theta, tau, phi_ny, theta_ny)
        position_vektor += np.array([dx, dy, dz])
        print('dz', dz)

        # Håll reda på ifall partikeln befinner sig i sfären eller inte.
        if np.sqrt(np.dot(position_vektor, position_vektor)) > radie_sfär:
            break

        # Plottvärderna för att se dosfördelningen, men får bara ut startpositionen inte alla andra delsteg...
        dos.append(energi - energi_efter_energiförlust(energi, steglängd, rho_medium, stopping_power_data))
        x.append(position_vektor[0])
        y.append(position_vektor[1])
        z.append(position_vektor[2])

        # Förlorar energi och får en ny efter endast steglängden
        energi = energi_efter_energiförlust(energi, steglängd, rho_medium, stopping_power_data)

        # Ändrar på vinklarna och värdet på steglängden innan nästa vinkeländring
        phi = phi_ny
        theta = theta_ny
        steglängd = tau

    # Beräkna den totala energideponeringen (i sfären) från partikeln.
    energideponering = energi_start - energi
    # print('Doslista', np.sum(dos))
    # print(f'energideponering: {energideponering} eV')
    return energideponering, x, y, z, dos
