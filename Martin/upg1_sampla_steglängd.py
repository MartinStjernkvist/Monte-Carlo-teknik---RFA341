from imports import *

@jit(nopython=True)
def invers_funktion(x, mu):
    """
    Invers funktion, som används för att sampla fotonens steglängd.
    :param x: Ett slumpat tal mellan 0 och 1.
    :param mu: Attenueringskoefficient
    """
    return -np.log(x) / mu

@jit(nopython=True)
def medelvägslängd(mu):
    """
    Funktion som samplar steglängden.
    Använder inverstransform-funktionen ovan.
    """

    # Steglängd angiven i voxelsidlängder.
    medelvägslängd = invers_funktion(np.random.rand(), mu) / voxel_sidlängd
    return medelvägslängd