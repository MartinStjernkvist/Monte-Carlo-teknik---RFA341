from imports import *

@jit(nopython=True)
def steg_vektoriserad(position_vektor, steglängd, phi, theta):

    # print(position_vektor)

    position_vektor = np.array([*position_vektor, 1], dtype=np.float64)

    steg_vektor = np.array([
        steglängd * np.sin(theta) * np.cos(phi),
        steglängd * np.sin(theta) * np.sin(phi),
        steglängd * np.cos(theta)
    ], dtype=np.float64)

    # print(steg_vektor)

    translation_matris = np.eye(4, dtype=np.float64)
    translation_matris[:3, 3] = steg_vektor

    ny_position_vektor = translation_matris @ position_vektor

    ny_position_vektor = ny_position_vektor[:3]

    # print(ny_position_vektor)


    return ny_position_vektor


@jit(nopython=True)
def steg(position_vektor, steglängd, phi, theta):
    """
    Funktion som tar ett steg i en specificerad riktning.
    :param theta: Spridningsvinkel.
    :param phi: Spridningsvinkel.
    :param x: Startposition x.
    :param y: Startposition y.
    :param z: Startposition z.
    :return: Ny position.
    """
    x = position_vektor[0]
    y = position_vektor[1]
    z = position_vektor[2]

    # Steg.
    dx = steglängd * np.sin(theta) * np.cos(phi) / voxel_sidlängd
    dy = steglängd * np.sin(theta) * np.sin(phi) / voxel_sidlängd
    dz = steglängd * np.cos(theta) / voxel_sidlängd

    # Tar steget.
    x = x + dx
    y = y + dy
    z = z + dz

    ny_position_vektor = np.array([x,y,z], dtype=np.float64)
    return ny_position_vektor


if __name__ == "__main__":
    position_vektor = np.array([1, 1, 1])
    steglängd = 5
    phi = pi
    theta = phi/4

    antal_iterationer = 10**3

    start = time.time()

    dummy_steg = steg_vektoriserad(position_vektor, steglängd, phi, theta)

    for i in range(antal_iterationer):
        steg = steg_vektoriserad(position_vektor, steglängd, phi, theta)
        # print(steg_vektoriserad(position_vektor, steglängd, phi, theta))

    end_time(start)

    start = time.time()

    dummy_steg_2 = steg(position_vektor, steglängd, phi, theta)

    for i in range(antal_iterationer):
        steg = steg(position_vektor, steglängd, phi, theta)
        # print(steg(position_vektor, steglängd, phi, theta))

    end_time(start)

