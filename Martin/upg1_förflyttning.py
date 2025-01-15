from imports import *

@jit(nopython=True)
def förflyttning(x, y, z, dx, dy, dz):
    x_ny = x + dx
    y_ny = y + dy
    z_ny = z + dz

    x_round = round(x_ny)
    y_round = round(y_ny)
    z_round = round(z_ny)

    return x_ny, y_ny, z_ny, x_round, y_round, z_round
