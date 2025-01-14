from imports import *


# @jit(nopython=True)
def ny_steg_transformera_koordinatsystem_3d(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B):
    """
    - Börjar på position A[x,y,z]- kalla detta koordinatsystem A.
    - Tar ett steg med steglängd steg_A_B, riktning (phi_A, theta_A), enligt koordinatsystemet i A.
            (exempelvis kommer phi_A att vara relativt enhetsvektorn i x-led för koord-syst A)
    - Tar ett steg till ny punkt - kalla denna punkt B.
    - Transformerar koordinatsystemet så att riktningsvektorn sammandfaller med
    nya koordinatsystemet.
            (nya enhetsvektorn i x-led, i B's koord-syst, ska ha samma riktning som fotonen
            hade när den tog steget)
    - Detta görs för att kunna sampla en ny riktning i nästa växelverkanprocess,
    då behövs nämligen ett koordinatsystem i B.

    :param steg_A_B: magnitud på steg från A till B
    :param phi_A: vinkel för steget mellan A och B
    :param theta_A: vinkel för steget mellan A och B

    :param steg_B_C: magnitud på steg från B till C
    :param phi_B: vinkel för steget mellan B och C
    :param theta_B: vinkel för steget mellan B och C

    :return: 3 värden är förflyttningen från B till C, enligt A's koord-syst
    """

    dx_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    dy_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    dz_A_B = steg_A_B * np.cos(theta_A)

    dx_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    dy_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    dz_B_C = steg_B_C * np.cos(theta_B)

    # transformera med rotation i z-led (phi)
    R_z = np.array(
        [
            [np.cos(phi_A), -np.sin(phi_A), 0],
            [np.sin(phi_A), np.cos(phi_A), 0],
            [0, 0, 1]
        ], dtype=np.float64)

    # för att x-axeln ska sammanfalla med riktningsvektorn måste rotationsvinkeln vara pi/2 - theta
    angle = -(pi / 2 - theta_A)
    R_y = np.array(
        [
            [np.cos(angle), 0, np.sin(angle)],
            [0, 1, 0],
            [-np.sin(angle), 0, np.cos(angle)],
        ], dtype=np.float64)

    # först rotation i theta (y-axeln), sedan rotation i phi (z-axeln)
    # R = R_z @ R_y
    R = np.dot(R_z, R_y)

    # Homogenous_matrix = np.array(
    #     [
    #         [R[0, 0], R[0, 1], R[0, 2], dx_A_B],
    #         [R[1, 0], R[1, 1], R[1, 2], dy_A_B],
    #         [R[2, 0], R[2, 1], R[2, 2], dz_A_B],
    #         [0, 0, 0, 1]
    #     ], dtype=np.float64)

    Homogenous_matrix = np.eye(4, dtype=np.float64)
    Homogenous_matrix[:3, :3] = R
    Homogenous_matrix[:3, 3] = np.array([dx_A_B, dy_A_B, dz_A_B], dtype=np.float64)

    # vektor_A_C = Homogenous_matrix @ np.array(
    #     [
    #         [dx_B_C],
    #         [dy_B_C],
    #         [dz_B_C],
    #         [1]   # nödvändigt för beräkningen
    #     ], dtype=np.float64)

    vektor_A_C = np.dot(Homogenous_matrix, np.array([dx_B_C, dy_B_C, dz_B_C, 1.0], dtype=np.float64))

    # vektor_A_B = np.array([
    #     [dx_A_B],
    #     [dy_A_B],
    #     [dz_A_B],
    #     [1]  # nödvändigt för beräkningen
    # ], dtype=np.float64)

    vektor_A_B = np.array([dx_A_B, dy_A_B, dz_A_B, 1.0], dtype=np.float64)

    # Vill ha vektor B->C
    vektor = vektor_A_C - vektor_A_B
    dx, dy, dz = vektor[0], vektor[1], vektor[2]

    return dx, dy, dz
