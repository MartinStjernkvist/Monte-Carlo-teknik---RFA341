from imports import *


def polar_vinkel(steglängd, scattering_power):
    """
    Funktion som ger polarvinkeln theta för elektronspridning.
    :param steglängd: Steglängden innan spridningen.
    :param scattering_power: Scattering power från annan funktion.
    :return: Spridningsvinkeln theta (sfäriska koordinater).
    """

    steglängd_cm = steglängd * 10 ** 2

    R = np.random.random()  # Slumpmässig tal mellan 0-1

    T = scattering_power

    theta_ny_polar = np.sqrt(-T * steglängd_cm * np.log(1 - R))
    theta_ny = theta_ny_polar

    # theta_ny_polar = random.choice([1, -1]) * theta_ny_polar
    # print('theta polar', theta_ny_polar)
    # theta_ny = np.pi / 2 - theta_ny_polar

    print(f'polarvinkel: {theta_ny:.3f}')
    return theta_ny
