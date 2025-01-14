from imports import *
from upg2_steg_transformation import ny_steg_transformera_koordinatsystem_3d
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd
from upg2_scattering_power import scattering_power_från_energi
from upg2_steglängd_från_energi import steglängd_från_energi
from upg2_elektron_polarvinkel import polar_vinkel

def laddad_partikel_väg_elektron(energi_start, position_start, phi, theta, radie_sfär, rho_medium,
                                 stopping_power_data, scatter_power_data, energiförlust_faktor):
    """
    Funktion som följer alfapartikeln allteftersom den växelverkar i ett medium.
    :param radie_sfär: Radien av sfären för fördelningen.
    :return: Energideponeringen innanför sfären.
    """
    x, y, z, dos = [], [], [], []

    # Startposition och startenergi.
    position_vektor = position_start
    energi = energi_start

    # Steglängd tas så att en viss energiförlust sker.
    steglängd = steglängd_från_energi(energi, rho_medium, stopping_power_data, energiförlust_faktor)

    # Den nya energin för elektronen.
    energi_ny = energi * energiförlust_faktor

    dx = steglängd * np.sin(theta) * np.cos(phi)
    dy = steglängd * np.sin(theta) * np.sin(phi)
    dz = steglängd * np.cos(theta)

    position_vektor += np.array([dx, dy, dz])
    print(f'positionsvektor: {position_vektor}')

    # dos.append(energi - energi_ny)
    # x.append(position_vektor[0])
    # y.append(position_vektor[1])
    # z.append(position_vektor[2])


    while energi > energi_start * 10**(-10):

        # Steglängd tas så att en viss energiförlust sker.
        steglängd_ny = steglängd_från_energi(energi, rho_medium, stopping_power_data, energiförlust_faktor)

        # Den nya energin för elektronen.
        energi_ny = energi * energiförlust_faktor

        # Scattering power för den nya energin, efter steget.
        scattering_power = scattering_power_från_energi(energi_ny, scatter_power_data, rho_medium)

        # Polarvinkel utifrån scattering power.
        theta_ny = polar_vinkel(steglängd, scattering_power)
        phi_ny = np.random.random() * 2 * pi

        # Ändrar på positionsvektor efter att transformations matrisen
        dx, dy, dz = ny_steg_transformera_koordinatsystem_3d(steglängd, phi, theta, steglängd_ny, phi_ny, theta_ny)
        position_vektor += np.array([dx, dy, dz])

        # Håll reda på ifall partikeln befinner sig i sfären eller inte.
        if np.sqrt(np.dot(position_vektor, position_vektor)) > radie_sfär:
            break

        else:
            # Plottvärderna för att se dosfördelningen, men får bara ut startpositionen inte alla andra delsteg...
            dos.append(energi - energi_ny)
            x.append(position_vektor[0])
            y.append(position_vektor[1])
            z.append(position_vektor[2])

            # Ändrar på vinklarna och värdet på steglängden innan nästa vinkeländring
            phi = phi_ny
            theta = theta_ny
            steglängd = steglängd_ny
            energi = energi_ny
            continue

    # Beräkna den totala energideponeringen (i sfären) från partikeln.
    energideponering = energi_start - energi
    # print(f'energideponering: {energideponering} eV')
    return energideponering, x, y, z, dos
