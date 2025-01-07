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
    Funktion som samplar steglängden utifrån den inverstransformerade funktionen ovan.
    """
    medelvägslängd = invers_funktion(np.random.rand(), mu) / voxel_sidlängd # LÄGG TILL VOXELLÄNGD
    return medelvägslängd