from imports import *
from upg1_attenueringsdata import attenueringsdata
from upg1_visualisera_bin_fil import fantom_matris


class Partikel:
    def __init__(self, start_position, start_energi, df_stopping_power):
        self.position = np.array(start_position)
        self.energi = start_energi
        self.df_stopping_power = df_stopping_power

    def riktning_alpha(self):
        theta = np.arccos(-1 + 2 * np.random.rand())
        phi = 2 * pi * np.random.rand()
        return theta, phi

    def steglängd_alpha(self):
        print('WIP')
        medelvägslängd = 100
        return medelvägslängd

    def förflyttning(self, steg_vektor):
        self.position += steg_vektor

    def energiförlust_alpha(self, steg):
        print('WIP')

        energiförlust = self.energi * 0.1
        self.energi -= energiförlust

        if self.energi <= 0:
            self.energi = 0

    def return_energi(self):
        energi = self.energi
        return energi

# @jit(nopython=True)
def laddad_partikel_väg_class(partikel, radie, max_antal_steg=100):

    start_energi = partikel.energi
    position_vektor = partikel.position
    steglängd = partikel.steglängd_alpha()
    theta, phi = partikel.riktning_uniform()

    print('start energi', start_energi)
    print('position_vektor', position_vektor)
    print('steglängd', steglängd)
    print('theta, phi', theta, phi)

    trajectory = [tuple(position_vektor)]

    steg_storlek = steglängd / max_antal_steg

    riktning = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.cos(phi), np.cos(theta)])
    riktning /= np.linalg.norm(riktning)
    steg_vektor = riktning * steg_storlek

    innanför = True
    summa_steglängd = 0

    sista_energi = partikel.return_energi()

    # Under tiden som partikeln fortfarnade inte tagit hela sitt steg.
    for i in range(max_antal_steg):

        # Medan partikeln fortfarande är innanför tumörcellen.
        while innanför == True:

            sista_energi = partikel.return_energi()  # Det sista värdet kommer sparas, resterande är irrelevanta.

            partikel.förflyttning(steg_vektor)
            partikel.energi_efter_energiförlust(steg_storlek)

            if np.dot(partikel.position, partikel.position) <= radie:
                innanför = True
                trajectory.append(tuple(partikel.position))
                print(f'Energideponering i position [{partikel.position}].')
            else:
                innanför = False
                print('Partikel utanför sfär!')

    energideponering = start_energi - sista_energi

    return trajectory, energideponering


if __name__ == "__main__":
    # Alpha partikel.
    df_stopping_power = pd.read_excel(attenueringsdata_file)

    partikel = Partikel((1.0, 1.0, 1.0), 10.0, df_stopping_power)

    radie = 4 # byt
    trajectory, energideponering = laddad_partikel_väg_class(partikel, radie, max_antal_steg=10000)

    # for step, pos in enumerate(trajectory):
    #     print(f"Step {step}: Position {pos}")
    print(f'\nTotal energideponering i sfär: {energideponering:.2f}.')

"""
Gammal fil: 

from imports import *
from upg1_attenueringsdata import attenueringsdata
from upg1_visualisera_bin_fil import fantom_matris

df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)
df_anatomidefinitioner = pd.read_excel(anatomidefinitioner_file, index_col=None)


class Partikel:
    def __init__(self, start_position, start_energi):
        self.position = np.array(start_position)
        self.energi = start_energi

    def förflyttning(self, steg_vektor):
        self.position += steg_vektor

    def energiförlust(self, steg, medium):
        energiförlust = 0.1 * self.energi  # byt
        self.energi -= energiförlust

        if self.energi <= 0:
            self.energi = 0


def laddad_partikel_väg(start_position, steglängd, phi, theta, matris, df_attenueringsdata, df_anatomidefinitioner,
                 max_antal_steg=100):
    position_vektor = partikel.position
    x, y, z = position_vektor[0], position_vektor[1], position_vektor[2]
    x_round, y_round, z_round = round(x), round(y), round(z)

    voxelvärde = matris[x_round, y_round, z_round]

    instans = attenueringsdata(voxelvärde, partikel.energi, df_attenueringsdata, df_anatomidefinitioner)
    start_medium = instans.voxelvärde_till_material()

    nuvarande_medium = start_medium

    trajectory = [tuple(partikel.position)]

    for _ in range(max_antal_steg):

        if partikel.energi <= 0:
            print(f'partikeln har slut på energi')
            break

        instans = attenueringsdata(voxelvärde, partikel.energi, df_attenueringsdata, df_anatomidefinitioner)
        nuvarande_medium = instans.voxelvärde_till_material()

        medelvägslängd = partikel.energi / (voxelvärde * 0.01)  # bara för att ha något
        steg_storlek = min(medelvägslängd, 1)

        riktning = np.random.normal(size=3)
        riktning /= np.linalg.norm(riktning)
        steg_vektor = riktning * steg_storlek

        partikel.förflyttning(steg_vektor)
        partikel.energiförlust(steg_storlek, voxelvärde)

        trajectory.append(tuple(partikel.position))

    return trajectory


if __name__ == "__main__":
    partikel = Partikel((1.0, 1.0, 1.0), 10.0)

    trajectory = laddad_partikel_väg(partikel, 100, fantom_matris)

    for step, pos in enumerate(trajectory):
        print(f"Step {step}: Position {pos}")
"""
