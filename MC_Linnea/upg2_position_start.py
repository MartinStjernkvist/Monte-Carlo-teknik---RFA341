from imports import *

@jit(nopython=True)
def position_start_innanför(radie_sfär):
    r = radie_sfär * np.random.rand()

    x = r
    # Utnyttja sfärisk symmetri.
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
def position_start_skal(radie_sfär, radie_partikel):
    r = radie_sfär - 0.5 * radie_partikel  # För att inte endast theta = pi ska ge utslag
    #Varför halva radien, räcker inte att ta hela radien??

    x = r
    # Utnyttja sfärisk symmetri.
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor