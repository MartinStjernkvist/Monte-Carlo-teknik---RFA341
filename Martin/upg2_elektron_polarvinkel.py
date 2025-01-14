from imports import *

def polar_vinkel(steglängd, scattering_power):

    steglängd_cm = steglängd * 10**2

    R = np.random.random()  # Slumpmässig tal mellan 0-1

    T = scattering_power

    theta_ny_polar = np.sqrt(-T * steglängd_cm * np.log(1 - R))
    theta_ny = theta_ny_polar
    print(f'polarvinkel: {theta_ny / (2 * np.pi)} * 2 pi')


    # theta_ny_polar = random.choice([1, -1]) * theta_ny_polar
    #
    # print('theta polar', theta_ny_polar)
    #
    # theta_ny = np.pi / 2 - theta_ny_polar
    #
    # print('theta ny', theta_ny)

    return theta_ny