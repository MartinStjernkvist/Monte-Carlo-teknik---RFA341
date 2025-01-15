from imports import *


def polar_vinkel(steglängd, scattering_power):
    """
    Funktion som ger polarvinkeln theta för elektronspridning.
    :param steglängd: Steglängden innan spridningen.
    :param scattering_power: Scattering power från annan funktion.
    :return: Spridningsvinkeln theta (sfäriska koordinater).
    """

    # steglängd_cm = steglängd * 10 ** 2

    R = np.random.random()  # Slumpmässig tal mellan 0-1

    T = scattering_power

    # Beräkning av polarvinkel, med inverstransform.
    theta_ny = np.sqrt(-T * steglängd * np.log(1 - R))

    return theta_ny
