from imports import *

def polar_vinkel(steglängd, scattering_power):

    R = np.random.random()  # Slumpmässig tal mellan 0-1

    T = scattering_power

    theta_ny_polar = np.sqrt(-T * steglängd * np.log(1 - R))

    theta_ny_polar = random.choice([1, -1]) * theta_ny_polar

    print('theta polar', theta_ny_polar)

    theta_ny_sfär = np.pi / 2 - theta_ny_polar

    print('theta ny', theta_ny_sfär)

    return theta_ny_sfär