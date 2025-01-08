from imports import *

@jit(nopython=True)
def förflyttning(x, y, z, dx, dy, dz):
    x = x + dx
    y = y + dy
    z = z + dz

    x_round = round(x)
    y_round = round(y)
    z_round = round(z)

    return x, y, z, x_round, y_round, z_round
