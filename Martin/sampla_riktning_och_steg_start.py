from imports import *

@jit(nopython=True)
def riktning_uniform():
    """
    Funktion som samplar uniformt spridningsvinklarna för en foton.
    Används både till startpositionen, samt vid fotoabsorption.
    :return: Spridningsvinklar theta och phi (sfäriska koordinater).
    """

    # cos(theta) ska vara mellan -1 och 1
    theta = np.arccos(-1 + 2 * random.rand())

    phi = 2 * pi * random.rand()
    return theta, phi

@jit(nopython=True)
def steg(theta, phi, steglängd, x, y, z):
    """
    Funktion som tar ett steg i en specificerad riktning.
    :param theta: Spridningsvinkel.
    :param phi: Spridningsvinkel.
    :param x: Startposition x.
    :param y: Startposition y.
    :param z: Startposition z.
    :return: Ny position.
    """

    # Steg.
    dx = steglängd * np.sin(theta) * np.cos(phi) / voxel_sidlängd
    dy = steglängd * np.sin(theta) * np.sin(phi) / voxel_sidlängd
    dz = steglängd * np.cos(theta) / voxel_sidlängd

    # Tar steget.
    x = x + dx
    y = y + dy
    z = z + dz

    return x, y, z