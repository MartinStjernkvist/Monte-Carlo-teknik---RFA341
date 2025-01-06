from ..imports import *
from ..attenueringsdata import attenueringsdata
from ..visualisera_bin_fil import fantom_matris

df_attenueringsdata = pd.read(attenueringsdata_file)
df_anatomidefinitioner = pd.read(anatomidefinitioner_file)


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


def laddad_partikel_väg(partikel, max_antal_steg, matris):
    position_vektor = partikel.position
    x, y, z = position_vektor[0], position_vektor[1], position_vektor[2]
    x_round, y_round, z_round = round(x), round(y), round(z)

    voxelvärde = matris[x_round, y_round, z_round]
    start_medium = attenueringsdata().voxelvärde_till_material(voxelvärde, partikel.energi, df_attenueringsdata,
                                                               df_anatomidefinitioner)
    nuvarande_medium = start_medium

    trajectory = [tuple(partikel.position)]

    for _ in range(max_antal_steg):

        if partikel.energi <= 0:
            print(f'partikeln har slut på energi')
            break

        nuvarande_medium = attenueringsdata().voxelvärde_till_material(voxelvärde, partikel.energi, df_attenueringsdata,
                                                                       df_anatomidefinitioner)

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
    partikel = Partikel((1, 1, 1), 10.0)

    trajectory = laddad_partikel_väg(partikel, fantom_matris)
