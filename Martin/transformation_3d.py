from imports import *


def transformera_koordinatsystem(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B):
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

    :return: 1) vektor vars första 3 värden är positionen för punkt C enligt A's koord-syst
            2) matris med enhetsvektorerna som B's koord-syst består av
    """

    dx_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    dy_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    dz_A_B = steg_A_B * np.cos(theta_A)

    dx_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    dy_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    dz_B_C = steg_B_C * np.cos(theta_B)

    enhets_vektorer_A = np.array(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 0]
        ])

    # transformera först med rotation i z-led (phi)
    # 4D matrix för att kunna skapa homogenitet-transformsmatris
    # R_z = np.array(
    #     [
    #         [np.cos(phi), -np.sin(phi), 0, 0],
    #         [np.sin(phi), np.cos(phi), 0, 0],
    #         [0, 0, 1, 0],
    #         [0, 0, 0, 0]
    #     ])

    R_z = np.array(
        [
            [np.cos(phi_A), -np.sin(phi_A), 0],
            [np.sin(phi_A), np.cos(phi_A), 0],
            [0, 0, 1]
        ])

    # transformera sedan i theta
    ### tänk att theta rotationsmatrisen har anti-clockwise rotation som standard, men theta är clockwise, så ta -theta
    # om vi har theta, måste kompenseringen (för att x-axeln ska hamna längs vektorn) vara pi/2 - theta, clockwise
    # R_y = np.array(
    #     [
    #         [np.cos(pi / 2 - theta), 0, -np.sin(pi / 2 - theta), 0],
    #         [0, 1, 0, 0],
    #         [np.sin(pi / 2 - theta), 0, np.cos(pi / 2 - theta), 0],
    #         [0, 0, 0, 0]
    #     ])

    angle = pi / 2 - theta_A
    R_y = np.array(
        [
            [np.cos(angle), 0, -np.sin(angle)],
            [0, 1, 0],
            [np.sin(angle), 0, np.cos(angle)],
        ])

    # R_y = np.identity(3)
    # R_z = np.identity(3)

    # R = R_y @ R_z

    # först rotation i theta (y-axeln), sedan rotation i phi (z-axeln)
    R = R_z @ R_y

    # print(f'\nR matrix: \n{R}')
    # print(f'\nlength of vector: \nx: {np.linalg.norm(R[0:3,0])}\ny: {np.linalg.norm(R[0:3,1])}\nz: {np.linalg.norm(R[0:3,2])}')

    Homogenous_matrix = np.array(
        [
            [R[0, 0], R[0, 1], R[0, 2], dx_A_B],
            [R[1, 0], R[1, 1], R[1, 2], dy_A_B],
            [R[2, 0], R[2, 1], R[2, 2], dz_A_B],
            [0, 0, 0, 1]
        ])

    vektor_A_C = Homogenous_matrix @ np.array(
        [
            [dx_B_C],
            [dy_B_C],
            [dz_B_C],
            [1]
        ])

    enhets_vektorer_B = Homogenous_matrix @ enhets_vektorer_A

    return vektor_A_C, enhets_vektorer_B


if __name__ == "__main__":

    #   ----------------------------------------------------------------------
    #   INPUT
    #   ---------------------------------------------------------------------

    # startpunkt A
    x_start = 1
    y_start = 1
    z_start = 1

    x_A = x_start
    y_A = y_start
    z_A = z_start

    # steg 1: från A till B
    theta_A = 3 * pi / 8
    phi_A = 1 * pi / 8
    steg_A_B = 3

    # steg 2: Från B till C, enligt koordinatsystemet för B
    theta_B = pi / 3
    phi_B = pi / 3
    steg_B_C = 2

    #   ----------------------------------------------------------------------
    #   BERÄKNING
    #   ---------------------------------------------------------------------

    dx_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    dy_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    dz_A_B = steg_A_B * np.cos(theta_A)

    dx_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    dy_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    dz_B_C = steg_B_C * np.cos(theta_B)

    vektor_A_C, enhets_vektorer_B = transformera_koordinatsystem(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B)

    x_A_C = vektor_A_C[0]
    x_tot = x_A + x_A_C

    y_A_C = vektor_A_C[1]
    y_tot = y_A + y_A_C

    z_A_C = vektor_A_C[2]
    z_tot = z_A + z_A_C

    #   ----------------------------------------------------------------------
    #   FELSÖKNING
    #   ---------------------------------------------------------------------

    print(f'tot: {np.array([[x_tot], [y_tot], [z_tot]])}')

    print(f'\nresult: \n{vektor_A_C}\nenhets_vektorer_B: \n{enhets_vektorer_B}')

    print(f'\ndot product x ({enhets_vektorer_B[0:3, 0]}), y ({enhets_vektorer_B[0:3, 1]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 1])}')
    print(f'dot product y ({enhets_vektorer_B[0:3, 1]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 1], enhets_vektorer_B[0:3, 2])}')
    print(f'dot product x ({enhets_vektorer_B[0:3, 0]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 2])}')



    #   ----------------------------------------------------------------------
    #   PLOTTNING
    #   ---------------------------------------------------------------------

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection='3d')
    ax.set_xlim([0, 6])
    ax.set_ylim([0, 6])
    ax.set_zlim([0, 6])

    ax.scatter(x_start, y_A, z_A, label='A', color='black', s=100)
    ax.scatter(x_start + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, label='B', color='grey', s=100)
    ax.scatter(x_A +x_A_C, y_A +y_A_C, z_A +z_A_C, label='C', color='purple', s=100)

    ax.quiver(x_A, y_A, z_A, 1, 0, 0, color='orange')
    ax.quiver(x_A, y_A, z_A, 0, 1, 0, color='brown')
    ax.quiver(x_A, y_A, z_A, 0, 0, 1, color='green')

    ax.quiver(x_A +dx_A_B, y_A +dy_A_B, z_A +dz_A_B, enhets_vektorer_B[0, 0], enhets_vektorer_B[1, 0], enhets_vektorer_B[2, 0],
              color='orange')
    ax.quiver(x_A +dx_A_B, y_A +dy_A_B, z_A +dz_A_B, enhets_vektorer_B[0, 1], enhets_vektorer_B[1, 1], enhets_vektorer_B[2, 1],
              color='brown')
    ax.quiver(x_A +dx_A_B, y_A +dy_A_B, z_A +dz_A_B, enhets_vektorer_B[0, 2], enhets_vektorer_B[1, 2], enhets_vektorer_B[2, 2],
              color='green')


    ax.quiver(x_A, y_A, z_A, dx_A_B, dy_A_B, dz_A_B, color='blue', label='första steget, A-> B')
    ax.quiver(x_A + dx_A_B, y_A +dy_A_B, z_A +dz_A_B, (x_A_C - dx_A_B), (y_A_C - dy_A_B), (z_A_C - dz_A_B), color='magenta',
              label='andra steget, B->C')
    ax.quiver(x_A, y_A, z_A, x_A_C, y_A_C, z_A_C, color='red')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.legend()
    plt.show()


    #   ----------------------------------------------------------------------
    #   INPUT - iteration 2
    #   ---------------------------------------------------------------------

    # startpunkt
    x_A = x_A + dx_A_B
    y_A = y_A + dy_A_B
    z_A = z_A + dz_A_B

    # steg 2
    steg_A_B = steg_B_C
    phi_A = phi_B
    theta_A = theta_B

    # steg 3
    steg_B_C = 2
    theta_B = pi / 4
    phi_B = pi / 4

    #   ----------------------------------------------------------------------
    #   BERÄKNING
    #   ---------------------------------------------------------------------

    dx_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    dy_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    dz_A_B = steg_A_B * np.cos(theta_A)

    dx_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    dy_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    dz_B_C = steg_B_C * np.cos(theta_B)

    vektor_A_C, enhets_vektorer_B = transformera_koordinatsystem(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B)

    x_A_C = vektor_A_C[0]
    x_tot = x_A + x_A_C

    y_A_C = vektor_A_C[1]
    y_tot = y_A + y_A_C

    z_A_C = vektor_A_C[2]
    z_tot = z_A + z_A_C

    #   ----------------------------------------------------------------------
    #   FELSÖKNING
    #   ---------------------------------------------------------------------

    print(f'tot: {np.array([[x_tot], [y_tot], [z_tot]])}')

    print(f'\nresult: \n{vektor_A_C}\nenhets_vektorer_B: \n{enhets_vektorer_B}')

    print(
        f'\ndot product x ({enhets_vektorer_B[0:3, 0]}), y ({enhets_vektorer_B[0:3, 1]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 1])}')
    print(
        f'dot product y ({enhets_vektorer_B[0:3, 1]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 1], enhets_vektorer_B[0:3, 2])}')
    print(
        f'dot product x ({enhets_vektorer_B[0:3, 0]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 2])}')

    #   ----------------------------------------------------------------------
    #   PLOTTNING
    #   ---------------------------------------------------------------------

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection='3d')
    ax.set_xlim([0, 6])
    ax.set_ylim([0, 6])
    ax.set_zlim([0, 6])

    ax.scatter(x_A, y_A, z_A, label='A', color='black', s=100)
    ax.scatter(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, label='B', color='grey', s=100)
    ax.scatter(x_A + x_A_C, y_A + y_A_C, z_A + z_A_C, label='C', color='purple', s=100)

    # ax.quiver(x_A, y_A, z_A, 1, 0, 0, color='orange')
    # ax.quiver(x_A, y_A, z_A, 0, 1, 0, color='brown')
    # ax.quiver(x_A, y_A, z_A, 0, 0, 1, color='green')

    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 0], enhets_vektorer_B[1, 0],
              enhets_vektorer_B[2, 0],
              color='orange')
    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 1], enhets_vektorer_B[1, 1],
              enhets_vektorer_B[2, 1],
              color='brown')
    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, enhets_vektorer_B[0, 2], enhets_vektorer_B[1, 2],
              enhets_vektorer_B[2, 2],
              color='green')

    ax.quiver(x_A, y_A, z_A, dx_A_B, dy_A_B, dz_A_B, color='blue', label='första steget, A-> B')
    ax.quiver(x_A + dx_A_B, y_A + dy_A_B, z_A + dz_A_B, (x_A_C - dx_A_B), (y_A_C - dy_A_B), (z_A_C - dz_A_B),
              color='magenta',
              label='andra steget, B->C')
    ax.quiver(x_A, y_A, z_A, x_A_C, y_A_C, z_A_C, color='red')
    ax.quiver(0, 0, 0, x_tot, y_tot, z_tot, color='black')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.legend()
    plt.show()
