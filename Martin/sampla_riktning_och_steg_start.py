from imports import *

@jit(nopython=True)
def riktning_koherent():
    # theta = random.gauss(pi / 2, 1)  # theta = np.arcsin(-1 + 2 * random.rand())

    theta = np.arccos(-1 + 2 * random.rand()) # cos(theta) ska vara mellan -1 och 1
    phi = 2 * pi * random.rand()
    return theta, phi

@jit(nopython=True)
def steg(theta, phi, steglängd, x, y, z):
    dx = steglängd * np.sin(theta) * np.cos(phi) / voxel_sidlängd
    dy = steglängd * np.sin(theta) * np.sin(phi) / voxel_sidlängd
    dz = steglängd * np.cos(theta) / voxel_sidlängd

    # new position
    x = x + dx
    y = y + dy
    z = z + dz
    return x, y, z