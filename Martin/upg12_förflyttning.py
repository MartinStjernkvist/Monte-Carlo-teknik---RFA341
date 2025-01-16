from imports import *


@jit(nopython=True)
def förflyttning(x, y, z, dx, dy, dz, voxel_sidlängd=1):
    x_ny = x + dx / voxel_sidlängd
    y_ny = y + dy / voxel_sidlängd
    z_ny = z + dz / voxel_sidlängd

    x_round = round(x_ny)
    y_round = round(y_ny)
    z_round = round(z_ny)

    return x_ny, y_ny, z_ny, x_round, y_round, z_round
