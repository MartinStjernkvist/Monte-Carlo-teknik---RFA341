from imports import *
from matris_njure import sliced_array_njure
from matris_benmärg import sliced_array_benmärg

x_size, y_size, z_size = sliced_array_njure.shape

deposited_energy = np.zeros((x_size, y_size, z_size))


def energideponering_benmärg(x_pos, y_pos, z_pos, energi):
    # antag att positionerna ges i term av voxelstorlek, alltså inte cm eller mm
    # -> vid stegberäkning måste stegen omvandlas till voxelstorlek innan de tas
    x_round = round(x_pos)
    y_round = round(y_pos)
    z_round = round(z_pos)

    if sliced_array_njure[x_round, y_round, z_round] != 0:
        return deposited_energy[x_round, y_round, z_round].append(energi)
