from imports import *


def riktning_start():
    # theta = random.gauss(pi / 2, 1)  # theta = np.arcsin(-1 + 2 * random.rand())

    theta = np.arccos(-1 + 2 * random.rand()) # cos(theta) ska vara mellan -1 och 1
    phi = 2 * pi * random.rand()
    return theta, phi
