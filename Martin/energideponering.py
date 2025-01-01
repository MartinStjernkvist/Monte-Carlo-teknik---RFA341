from imports import *
from matriser import slicad_njure_matris, slicad_benmärg_matris

x_size, y_size, z_size = slicad_njure_matris.shape

benmärg_matris_deponerad_energi = np.zeros((x_size, y_size, z_size))


def energideponering_benmärg(x_pos, y_pos, z_pos, energi, slicad_benmärg_matris):
    # antag att positionerna ges i term av voxelstorlek, alltså inte cm eller mm
    # -> vid stegberäkning måste stegen omvandlas till voxelstorlek innan de tas
    x_round = round(x_pos)
    y_round = round(y_pos)
    z_round = round(z_pos)

    if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
        benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energi
    return benmärg_matris_deponerad_energi


if __name__ == "__main__":
    start = time.time()

    test_matris = np.zeros((5, 5))
    print(test_matris)
    test_matris[0,0] += 1
    print(test_matris)
    test_matris[0,0] += 1
    print(test_matris)

    end_time(start)