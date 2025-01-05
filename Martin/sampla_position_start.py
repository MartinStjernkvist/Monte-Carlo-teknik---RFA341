from imports import *

@jit(nopython=True)
def position_start(slicad_njure_matris):
    """
    Funktion som samplar en startposition utifrån en matris med endast njurvävnad.
    :param slicad_njure_matris: Matris med njurvävnad.
    :return: Slumpad startposition
    """
    x_size, y_size, z_size = slicad_njure_matris.shape

    # Tar första möjliga voxel vars värde inte är 0.
    while True:
        x = np.random.randint(0, x_size)
        y = np.random.randint(0, y_size)
        z = np.random.randint(0, z_size)

        if slicad_njure_matris[x, y, z] != 0:
            break

    return x, y, z

if __name__ == "__main__":
    from matriser import slicad_njure_matris

    x,y,z = position_start(slicad_njure_matris)
    print(x,y,z)