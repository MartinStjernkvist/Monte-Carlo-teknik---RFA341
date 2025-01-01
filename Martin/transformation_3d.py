from imports import *


def transform(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B):
    """
    - Börjar på x,y,z - kalla detta koordinatsystem A.
    - Tar ett steg med steglängd s, riktning (phi, theta), enligt koordinatsystemet i A.
            (exempelvis kommer phi att vara relativt enhetsvektorn i x-led för koord-syst A)
    - Tar ett steg till ny punkt - kalla denna punkt B
    - Transformerar koordinatsystemet så att riktningsvektorn sammandfaller med
    nya koordinatsystemet
            (nya enhetsvektorn i x-led, i B's koord-syst, ska ha samma riktning som fotonen
            hade när den tog steget)
    - Detta görs för att kunna sampla en ny riktning i nästa växelverkanprocess,
    då behövs nämligen ett koordinatsystem i B som
    :param steg_A_B: magnitud på steg från A till B
    :param phi_A: vinkel för steget mellan A och B
    :param theta_A: vinkel för steget mellan A och B
    :param steg_B_C: magnitud på steg från B till C
    :param phi_B: vinkel för steget mellan B och C
    :param theta_B: vinkel för steget mellan B och C
    :return: 1) resultatvektor, vars första 3 värden är positionen för punkt C enligt A's koord-syst
            2) matris med enhetsvektorerna som B's koord-syst består av
    """

    # """
    # - Börjar på x,y,z - kalla detta koordinatsystem A.
    # - Tar ett steg med steglängd s, riktning (phi, theta), enligt koordinatsystemet i A.
    #         (exempelvis kommer phi att vara relativt enhetsvektorn i x-led för koord-syst A)
    # - Tar ett steg till ny punkt - kalla denna punkt B
    # - Transformerar koordinatsystemet så att riktningsvektorn sammandfaller med
    # nya koordinatsystemet
    #         (nya enhetsvektorn i x-led, i B's koord-syst, ska ha samma riktning som fotonen
    #         hade när den tog steget)
    # - Detta görs för att kunna sampla en ny riktning i nästa växelverkanprocess,
    # då behövs nämligen ett koordinatsystem i B som
    # :param x_B: x-komponent av ny riktningsvektor (steg) i B's koordinatsystem
    # :param y_B: y-komponent --''--
    # :param z_B: z-komponent --''--
    # :param d_x: x-komp av steg från A till B
    # :param d_y: y-komp --''--
    # :param d_z: z-komp --''--
    # :param phi_A: vinkel sfärisk
    # :param theta_A: vinkel sfärisk
    # """

    d_x_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    d_y_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    d_z_A_B = steg_A_B * np.cos(theta_A)

    d_x_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    d_y_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    d_z_B_C = steg_B_C * np.cos(theta_B)


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
            [R[0, 0], R[0, 1], R[0, 2], d_x_A_B],
            [R[1, 0], R[1, 1], R[1, 2], d_y_A_B],
            [R[2, 0], R[2, 1], R[2, 2], d_z_A_B],
            [0, 0, 0, 1]
        ])

    result = Homogenous_matrix @ np.array(
        [
            [d_x_B_C],
            [d_y_B_C],
            [d_z_B_C],
            [1]
        ])

    enhets_vektorer_B = Homogenous_matrix @ enhets_vektorer_A
    return result, enhets_vektorer_B


if __name__ == "__main__":
    start = time.time()

    # steg 1: från A till B
    theta_A = 3 * pi / 8
    phi_A = 3 * pi / 8
    steg_A_B = 3
    # d_x = steg * np.sin(theta) * np.cos(phi)
    # d_y = steg * np.sin(theta) * np.sin(phi)
    # d_z = steg * np.cos(theta)


    # steg 2: Från B till C, enligt koordinatsystemet för B
    theta_B = pi / 4
    phi_B = pi / 4
    steg_B_C = 1
    # x_B = 1
    # y_B = 1
    # z_B = 1

    d_x_A_B = steg_A_B * np.sin(theta_A) * np.cos(phi_A)
    d_y_A_B = steg_A_B * np.sin(theta_A) * np.sin(phi_A)
    d_z_A_B = steg_A_B * np.cos(theta_A)

    d_x_B_C = steg_B_C * np.sin(theta_B) * np.cos(phi_B)
    d_y_B_C = steg_B_C * np.sin(theta_B) * np.sin(phi_B)
    d_z_B_C = steg_B_C * np.cos(theta_B)

    result, enhets_vektorer_B = transform(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B)
    x_A = result[0]
    y_A = result[1]
    z_A = result[2]

    print(f'\nresult: \n{result}\nenhets_vektorer_B: \n{enhets_vektorer_B}')

    print(f'\ndot product x ({enhets_vektorer_B[0:3, 0]}), y ({enhets_vektorer_B[0:3, 1]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 1])}')
    print(f'dot product y ({enhets_vektorer_B[0:3, 1]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 1], enhets_vektorer_B[0:3, 2])}')
    print(f'dot product x ({enhets_vektorer_B[0:3, 0]}), z ({enhets_vektorer_B[0:3, 2]}): {np.dot(enhets_vektorer_B[0:3, 0], enhets_vektorer_B[0:3, 2])}')

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection='3d')
    ax.set_xlim([0, 4])
    ax.set_ylim([0, 4])
    ax.set_zlim([0, 4])

    ax.scatter(0, 0, 0, label='A', color='black', s=100)
    ax.scatter(d_x_A_B, d_y_A_B, d_z_A_B, label='B', color='grey', s=100)
    ax.scatter(x_A, y_A, z_A, label='C', color='purple', s=100)

    ax.quiver(0, 0, 0, 1, 0, 0, color='orange')
    ax.quiver(0, 0, 0, 0, 1, 0, color='brown')
    ax.quiver(0, 0, 0, 0, 0, 1, color='green')

    ax.quiver(d_x_A_B, d_y_A_B, d_z_A_B, enhets_vektorer_B[0, 0], enhets_vektorer_B[1, 0], enhets_vektorer_B[2, 0], color='orange')
    ax.quiver(d_x_A_B, d_y_A_B, d_z_A_B, enhets_vektorer_B[0, 1], enhets_vektorer_B[1, 1], enhets_vektorer_B[2, 1], color='brown')
    ax.quiver(d_x_A_B, d_y_A_B, d_z_A_B, enhets_vektorer_B[0, 2], enhets_vektorer_B[1, 2], enhets_vektorer_B[2, 2], color='green')

    ax.quiver(0, 0, 0, x_A, y_A, z_A, color='red')
    ax.quiver(0, 0, 0, d_x_A_B, d_y_A_B, d_z_A_B, color='blue')
    ax.quiver(d_x_A_B, d_y_A_B, d_z_A_B, (x_A - d_x_A_B), (y_A - d_y_A_B), (z_A - d_z_A_B), color='magenta')

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.legend()
    plt.show()

    end_time(start)
