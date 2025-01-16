from imports import *
from upg12_steg_transformation import transformera_koordinatsystem
from upg2_stopping_power_och_steglängd import stopping_power_och_steglängd
from upg2_scattering_power import scattering_power_från_energi
from upg2_steglängd_från_energi import steglängd_från_energi
from upg2_elektron_polarvinkel import polar_vinkel
from upg1_förflyttning import förflyttning
from upg12_rotation_matris import rotations_matris


def laddad_partikel_väg_elektron(energi_start, position_start, phi_start, theta_start, radie_sfär, rho_medium,
                                 stopping_power_data, scatter_power_data, energiförlust_faktor):
    """
    Funktion som följer alfapartikeln allteftersom den växelverkar i ett medium.
    :param energi_start: Elektronens startenergi.
    :param position_start: Elektronens startposition.
    :param phi_start: Sfärisk vinkel.
    :param theta_start: Sfärisk vinkel.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param rho_medium: Mediumets densitet (kg / m^3).
    :param stopping_power_data: Stopping power tabelldata.
    :param scatter_power_data: Scatter power tabelldata.
    :param energiförlust_faktor: Energiförlustfaktor efter varje steglängd.
    :return: Energideponeringen innanför sfären.
    """

    print(
        f'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nny partikel:\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    # Initiera tomma listor för att spara datan.
    x_list, y_list, z_list, dos = [], [], [], []

    # Startposition och startenergi.
    position_vektor = position_start
    x_start, y_start, z_start = position_vektor[0], position_vektor[1], position_vektor[2]

    energi = energi_start
    print(f'start energi: {energi:.3f} eV')

    # Steglängd tas så att en viss energiförlust sker.
    steglängd = steglängd_från_energi(energi, rho_medium, stopping_power_data, energiförlust_faktor)
    print(f'steglängd: {steglängd * 10 ** 6:.2f} mikrometer')

    # Den nya energin för elektronen.
    energi_ny = energi * energiförlust_faktor

    # Första steget.
    dx = steglängd * np.sin(theta_start) * np.cos(phi_start)
    dy = steglängd * np.sin(theta_start) * np.sin(phi_start)
    dz = steglängd * np.cos(theta_start)

    x, y, z, _, _, _ = förflyttning(x_start, y_start, z_start, dx, dy, dz)
    position_vektor = np.array([x, y, z])
    # print(f'positionsvektor: {position_vektor}')

    if np.sqrt(np.dot(position_vektor, position_vektor)) < radie_sfär:

        # Spara mätpunkter.
        dos.append(energi - energi_ny)
        x_list.append(x)
        y_list.append(y)
        z_list.append(z)

        # Ändrar på vinklarna, värdet på steglängden och energin.
        phi = phi_start
        theta = theta_start
        energi = energi_ny
        R = rotations_matris(phi, theta)

        print(
            f'\n-----------------------------------\ninitiera while loop:\n-----------------------------------')

        #   -----------------------------------
        #   Medan elektronen befinner sig i sfären
        #   och
        #   energin är över en tröskelenergi.
        #   -----------------------------------
        while np.sqrt(np.dot(position_vektor, position_vektor)) < radie_sfär and energi > energi_start * 10 ** (-6):
            print('')

            # Steglängd tas så att en viss energiförlust sker.
            steglängd_ny = steglängd_från_energi(energi, rho_medium, stopping_power_data, energiförlust_faktor)
            print(f'steglängd: {steglängd_ny * 10 ** 6:.2f} mikrometer')

            # Den nya energin för elektronen.
            energi_ny = energi * energiförlust_faktor
            print(f'energi ny: {energi:.2f} eV')

            # Scattering power för den nya energin, efter steget.
            scattering_power = scattering_power_från_energi(energi_ny, scatter_power_data, rho_medium)
            print(f'scattering power T: {scattering_power / 100:.2f} rad^2 / cm')

            # Polarvinkel utifrån scattering power.
            theta_ny = polar_vinkel(steglängd, scattering_power)
            print(f'polarvinkel: {theta_ny * 360 / (2 * np.pi):.2f} grader')

            phi_ny = np.random.random() * 2 * pi

            # Koordinattransformation ger ny positionsvektor.
            dx, dy, dz, R_ny = transformera_koordinatsystem(steglängd, phi, theta, steglängd_ny, phi_ny, theta_ny, R)

            x, y, z, _, _, _ = förflyttning(x, y, z, dx, dy, dz)
            position_vektor = np.array([x,y,z])

            #   -----------------------------------
            #   Håll reda på ifall:
            #   partikeln befinner sig i sfären.
            #   -----------------------------------
            if np.sqrt(np.dot(position_vektor, position_vektor)) > radie_sfär:
                print('\n!!!UTANFÖR!!!')
                break

            # Om elektronen fortfarande är i sfären -> spara mätpunkter.
            else:
                # Spara mätpunkter.
                dos.append(energi - energi_ny)
                x_list.append(x)
                y_list.append(y)
                z_list.append(z)

                # Ändrar på vinklarna, värdet på steglängden och energin.
                phi = phi_ny
                theta = theta_ny
                steglängd = steglängd_ny
                energi = energi_ny
                R = R_ny

                print(f'energi: {energi:.2f} eV')
                print('Rotationsmatris:\n', R)
                continue

    else:
        print('\n!!!UTANFÖR!!!')

    # Beräkna den totala energideponeringen (i sfären) från partikeln.
    energideponering = energi_start - energi
    # print(f'energideponering: {energideponering} eV')
    return energideponering, x_list, y_list, z_list, dos
