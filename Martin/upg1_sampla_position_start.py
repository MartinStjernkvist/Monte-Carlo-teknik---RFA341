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

        # x = np.random.rand(0, x_size)
        # y = np.random.rand(0, y_size)
        # z = np.random.rand(0, z_size)
        #
        # x_round = float(round(x))
        # y_round = float(round(y))
        # z_round = float(round(z))
        #
        # if slicad_njure_matris[x_round, y_round, z_round] != 0:
        #     break

    return x, y, z  # , x_round, y_round, z_round


if __name__ == "__main__":
    from upg1_matriser import slicad_njure_matris

    x, y, z = position_start(slicad_njure_matris)
    print(x, y, z)
