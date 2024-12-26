from imports import *

def steg(theta, phi, steglängd, x,y,z):
    voxel_sidlängd = 0.15 # cm
    dx = steglängd * np.cos(theta) * np.cos(phi) / voxel_sidlängd
    dy = steglängd * np.cos(theta) * np.sin(phi) / voxel_sidlängd
    dz = steglängd * np.sin(theta) / voxel_sidlängd

    # new position
    x = x + dx
    y = y + dy
    z = z + dz
    return x, y, z