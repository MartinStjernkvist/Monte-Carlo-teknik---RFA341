from imports import *

@jit(nopython=True)
def riktning_uniform():
    """
    Uniform riktning.
    """
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi


@jit(nopython=True)
def riktning_skal():
    """
    Ifall skalförfördelnning:
    pi / 2 < phi < 3 * pi / 2 för att effektivisera koden.
    Då kan antalet iterationer halveras, eftersom man vet att
    hälften av partiklarna ändå skulle lämna sfären.
    """
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = pi / 2 * (2*np.random.rand() + 1)
    return theta, phi