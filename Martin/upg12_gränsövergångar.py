from imports import *
from upg1_attenueringsdata import attenueringsdata
from upg1_visualisera_bin_fil import fantom_matris


def partikel_väg(start_position, steglängd, phi, theta, matris, df_attenueringsdata, df_anatomidefinitioner,
                 max_antal_steg=100):
    position_vektor = start_position
    # x, y, z = position_vektor[0], position_vektor[1], position_vektor[2]
    # x_round, y_round, z_round = round(x), round(y), round(z)

    trajectory = [tuple(start_position)]

    x_round = round(position_vektor[0])
    y_round = round(position_vektor[1])
    z_round = round(position_vektor[2])

    voxelvärde = matris[x_round, y_round, z_round]
    instans = attenueringsdata(voxelvärde, 0, df_attenueringsdata, df_anatomidefinitioner)
    start_material, _ = instans.voxelvärde_till_material()

    material_set = {start_material}

    steg_storlek = steglängd / max_antal_steg

    riktning = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.cos(phi), np.cos(theta)])
    riktning /= np.linalg.norm(riktning)
    steg_vektor = riktning * steg_storlek

    summa_steglängd = 0

    while abs(summa_steglängd - steglängd) > 10 ** -5:
        summa_steglängd += steg_storlek

        position_vektor += steg_vektor

        trajectory.append(tuple(position_vektor))

        x_round = round(position_vektor[0])
        y_round = round(position_vektor[1])
        z_round = round(position_vektor[2])

        voxelvärde = matris[x_round, y_round, z_round]
        instans = attenueringsdata(voxelvärde, 0, df_attenueringsdata, df_anatomidefinitioner)
        material, _ = instans.voxelvärde_till_material()

        material_set.add(material)

    return trajectory, material_set


if __name__ == "__main__":

    start_position = np.array([1, 1, 1], dtype=np.float64)

    steglängd = 400
    phi = pi / 4
    theta = pi / 5
    matris = fantom_matris
    df_attenueringsdata = pd.read_excel(attenueringsdata_file)
    df_anatomidefinitioner = pd.read_excel(anatomidefinitioner_file)

    trajectory, material_set = partikel_väg(start_position, steglängd, phi, theta, matris, df_attenueringsdata,
                              df_anatomidefinitioner, max_antal_steg=100)

    for step, pos in enumerate(trajectory):
        print(f"Step {step}: Position {pos}")

    print(material_set)
