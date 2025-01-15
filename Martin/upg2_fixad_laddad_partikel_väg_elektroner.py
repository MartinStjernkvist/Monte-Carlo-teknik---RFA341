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
    :param energi_start: Elektronens startenergi.
    :param position_start: Elektronens startposition.
    :param phi: Sfärisk vinkel.
    :param theta: Sfärisk vinkel.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param rho_medium: Mediumets densitet (kg / m^3).
    :param stopping_power_data: Stopping power tabelldata.
    :param scatter_power_data: Scatter power tabelldata.
    :param energiförlust_faktor: Energiförlustfaktor efter varje steglängd.
    :return: Energideponeringen innanför sfären.
    """

    print(f'\nny partikel:\n')

    # Initiera tomma listor för att spara datan.
    x, y, z, dos = [], [], [], []

    # Startposition och startenergi.
    position_vektor = position_start
    energi = energi_start
    print(f'start energi: {energi:.3f} eV')

    # Steglängd tas så att en viss energiförlust sker.
    steglängd = steglängd_från_energi(energi, rho_medium, stopping_power_data, energiförlust_faktor)
    print(f'steglängd: {steglängd * 10**6:.2f} mikrometer')


    # Den nya energin för elektronen.
    energi_ny = energi * energiförlust_faktor

    # Första steget.
    dx = steglängd * np.sin(theta) * np.cos(phi)
    dy = steglängd * np.sin(theta) * np.sin(phi)
    dz = steglängd * np.cos(theta)

    position_vektor += np.array([dx, dy, dz])
    # print(f'positionsvektor: {position_vektor}')

    if np.sqrt(np.dot(position_vektor, position_vektor)) < radie_sfär:
        dos.append(energi - energi_ny)
        x.append(position_vektor[0])
        y.append(position_vektor[1])
        z.append(position_vektor[2])

        print(f'\ninitiera while loop:')
        energi = energi_ny

    #   ----------------------------------------------------------------------
    #   Medan elektronen befinner sig i sfären
    #   och
    #   har en energi som är över en tröskelenergi.
    #   ----------------------------------------------------------------------

    while np.sqrt(np.dot(position_vektor, position_vektor)) < radie_sfär and energi > energi_start * 10 ** (-6):
        print('\n')

        # Steglängd tas så att en viss energiförlust sker.
        steglängd_ny = steglängd_från_energi(energi, rho_medium, stopping_power_data, energiförlust_faktor)
        print(f'steglängd: {steglängd_ny * 10 ** 6:.3f} mikrometer')

        # Den nya energin för elektronen.
        energi_ny = energi * energiförlust_faktor
        print(f'energi ny: {energi:.3f} eV')

        # Scattering power för den nya energin, efter steget.
        scattering_power = scattering_power_från_energi(energi_ny, scatter_power_data, rho_medium)
        print(f'scattering power T: {scattering_power / 100:.3f} rad^2 / cm')

        # Polarvinkel utifrån scattering power.
        theta_ny = polar_vinkel(steglängd, scattering_power)
        print(f'polarvinkel: {theta_ny * 360 / (2 * np.pi):.3f} grader')

        phi_ny = np.random.random() * 2 * pi

        # Ändrar på positionsvektor efter att transformations matrisen
        dx, dy, dz = ny_steg_transformera_koordinatsystem_3d(steglängd, phi, theta, steglängd_ny, phi_ny, theta_ny)
        position_vektor += np.array([dx, dy, dz])

        #   ----------------------------------------------------------------------
        #   Håll reda på ifall partikeln befinner sig i sfären eller inte.
        #   ----------------------------------------------------------------------
        if np.sqrt(np.dot(position_vektor, position_vektor)) > radie_sfär:
            break

        # Om elektronen fortfarande är i sfären -> spara mätpunkter.
        else:
            # Mätpunkter för att plotta dosfördelningen.
            dos.append(energi - energi_ny)
            x.append(position_vektor[0])
            y.append(position_vektor[1])
            z.append(position_vektor[2])

            # Ändrar på vinklarna, värdet på steglängden och energin innan nästa vinkeländring.
            phi = phi_ny
            theta = theta_ny
            steglängd = steglängd_ny
            energi = energi_ny

            print(f'energi: {energi:.3f} eV')
            continue

    # Beräkna den totala energideponeringen (i sfären) från partikeln.
    energideponering = energi_start - energi
    # print(f'energideponering: {energideponering} eV')
    return energideponering, x, y, z, dos
