from imports import *


def transformera_koordinatsystem(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B, R_A):
    """
    Koordinattransformation -> ny positionsvektor.

    1)  Börjar på position A[x,y,z]- kalla detta koordinatsystem A.
    2)  Tar ett steg med steglängd steg_A_B, riktning (phi_A, theta_A), enligt koordinatsystemet i A.
                - Vinklarna (phi_A, theta_A) är relativa koordinatsystemet i A.
    3)  Efter steget befinner sig partikeln i en ny punkt - kalla denna punkt B.
    4)  Transformerar koordinatsystemet så att riktningsvektorn sammanfaller med nya koordinatsystemet.
                - Nya enhetsvektorn i z-led, i B's koordinatsystem,
                ska ha samma riktning som fotonen hade när den tog steget
    5)  Detta görs för att kunna sampla en ny riktning i nästa växelverkanprocess,
        då behövs nämligen ett koordinatsystem i B.

    :param steg_A_B: Magnitud på steg från A till B.
    :param phi_A: Vinkel för steget mellan A och B.
    :param theta_A: Vinkel för steget mellan A och B.
    :param steg_B_C: Magnitud på steg från B till C.
    :param phi_B: Vinkel för steget mellan B och C.
    :param theta_B: Vinkel för steget mellan B och C.
    :return: Vektorkomponenter för vektor_BC (A).
    """

    #   -----------------------------------
    #   Vektorkomponenter.
    #   -----------------------------------
    x_B_A = steg_A_B * np.sin(theta_A) * np.cos(phi_A)  # x_B (A)
    y_B_A = steg_A_B * np.sin(theta_A) * np.sin(phi_A)  # y_B (A)
    z_B_A = steg_A_B * np.cos(theta_A)  # z_B (A)

    x_C_B = steg_B_C * np.sin(theta_B) * np.cos(phi_B)  # x_C (B)
    y_C_B = steg_B_C * np.sin(theta_B) * np.sin(phi_B)  # y_C (B)
    z_C_B = steg_B_C * np.cos(theta_B)  # z_C (B)

    #   -----------------------------------
    #   Vektorer.
    #   -----------------------------------

    vektor_B_A_innan = np.array([x_B_A, y_B_A, z_B_A])  # vektor_B (A)

    vektor_C_B = np.array([x_C_B, y_C_B, z_C_B])  # vektor_C (B)

    #   -----------------------------------
    #   Ifall en partikel åkt till A, måste ett virtuellt koordinatsystemet
    #   för A sammanfalla med riktningsvektorn för samma partikel.
    #
    #   Observera att koordinatsystemet för A egentligen kommer vara
    #   detsamma som det generella (t.ex. koordinatsystemet för fantomen
    #   eller sfären), men att en virtuell transformation måste
    #   göras för att kunna utföra beräkningar med en homogen matris.
    #   -----------------------------------

    # Rotera vektor_B_A enligt riktningen för föregående partikel.
    # Samma som att ett virtuellt koordinatsystem för A används.
    vektor_B_A = np.dot(R_A, vektor_B_A_innan)

    #   -----------------------------------
    #   Rotationsmatrisen: R.
    #   -----------------------------------

    # Rotation i z-led (phi).
    R_z = np.array(
        [
            [np.cos(phi_A), -np.sin(phi_A), 0],
            [np.sin(phi_A), np.cos(phi_A), 0],
            [0, 0, 1]
        ], dtype=np.float64)

    # Rotation i y-led (theta).
    # För att z-axeln ska sammanfalla med riktningsvektorn
    # -> måste rotationsvinkeln vara theta_A.
    R_y = np.array(
        [
            [np.cos(theta_A), 0, np.sin(theta_A)],
            [0, 1, 0],
            [-np.sin(theta_A), 0, np.cos(theta_A)],
        ], dtype=np.float64)

    # Först rotation i y-led, sedan rotation i z-led.
    R_innan = np.dot(R_z, R_y)

    # Rotera koordinatsystemet enligt föregående riktningsvektor.
    R = np.dot(R_A, R_innan)

    #   -----------------------------------
    #   Homogen matris: H.
    #   -----------------------------------
    H = np.eye(4, dtype=np.float64)
    H[:3, :3] = R
    H[:3, 3] = np.array([vektor_B_A[0], vektor_B_A[1], vektor_B_A[2]], dtype=np.float64)

    #   -----------------------------------
    #   Vektorer.
    #   -----------------------------------

    # vektor_C (A)
    vektor_C_A = np.dot(H, np.array([vektor_C_B[0], vektor_C_B[1], vektor_C_B[2], 1.0], dtype=np.float64))

    # vektor_B_A fast med en etta i slutet.
    vektor_B_A = np.array([vektor_B_A[0], vektor_B_A[1], vektor_B_A[2], 1.0])

    #   -----------------------------------
    #   Förflyttning från B -> C, i A's koordinatsystem.
    #   -----------------------------------

    # Vill ha vektor_BC (A).
    # Vektorn från B till C, i A's koordinatsystem.
    vektor_BC_A = vektor_C_A - vektor_B_A

    # Vektorkomponenter av vektor_BC (A).
    x_BC_A, y_BC_A, z_BC_A = vektor_BC_A[0], vektor_BC_A[1], vektor_BC_A[2]

    # Genom att dela upp i vektorkomponenter kan positionen
    # för partikeln följas i det övegripande koordinatsystemet,
    # genom att addera stegvektorer.
    dx, dy, dz = x_BC_A, y_BC_A, z_BC_A

    return dx, dy, dz, R
