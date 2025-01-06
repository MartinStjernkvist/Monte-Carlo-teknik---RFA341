

def laddad_partikel_steg(x,y,z, max_antal_steg, matris):
    x_round, y_round, z_round = round(x), round(y), round(z)
    start_medium = voxelvärde_till_material(x_round, y_round, z_round, matris)
    nuvarande_medium = start_medium

    while nuvarande_medium == start_medium:

        for i in range(max_antal_steg):

