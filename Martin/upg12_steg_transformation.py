from imports import *


def transformera_koordinatsystem(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B):
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
    :return: Vektorkomponenter för vektor_C.
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
    #   Rotationsmatriser: R.
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
    vinkel = theta_A
    R_y = np.array(
        [
            [np.cos(vinkel), 0, np.sin(vinkel)],
            [0, 1, 0],
            [-np.sin(vinkel), 0, np.cos(vinkel)],
        ], dtype=np.float64)

    # Först rotation i y-led, sedan rotation i z-led.
    R = np.dot(R_z, R_y)

    #   -----------------------------------
    #   Homogen matris: H.
    #   -----------------------------------
    H = np.eye(4, dtype=np.float64)
    H[:3, :3] = R
    H[:3, 3] = np.array([x_B_A, y_B_A, z_B_A], dtype=np.float64)

    #   -----------------------------------
    #   Vektorer.
    #   -----------------------------------

    # vektor_C (A)
    vektor_C_A = np.dot(H, np.array([x_C_B, y_C_B, z_C_B, 1.0], dtype=np.float64))

    # vektor_B (A)
    vektor_B_A = np.array([x_B_A, y_B_A, z_B_A, 1.0], dtype=np.float64)

    # Vill ha vektor_C.
    vektor = vektor_C_A - vektor_B_A

    # Vektorkomponenter av vektor_C.
    x_C, y_C, z_C = vektor[0], vektor[1], vektor[2]

    return x_C, y_C, z_C
