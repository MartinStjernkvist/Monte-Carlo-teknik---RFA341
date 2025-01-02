from imports import *


def position_start(slicad_njure_matris):
    x_size, y_size, z_size = slicad_njure_matris.shape

    while True:
        x = np.random.randint(0, x_size)
        y = np.random.randint(0, y_size)
        z = np.random.randint(0, z_size)

        if slicad_njure_matris[x, y, z] != 0:
            break

    return x, y, z
