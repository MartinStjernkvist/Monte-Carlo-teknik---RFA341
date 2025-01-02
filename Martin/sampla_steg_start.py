from imports import *


def första_steg(theta, phi, steglängd, x, y, z):
    dx = steglängd * np.cos(theta) * np.cos(phi) / voxel_sidlängd
    dy = steglängd * np.cos(theta) * np.sin(phi) / voxel_sidlängd
    dz = steglängd * np.sin(theta) / voxel_sidlängd

    # new position
    x = x + dx
    y = y + dy
    z = z + dz
    return x, y, z
