from imports import *


def rotations_matris(phi, theta):
    # Rotation i z-led (phi).
    vinkel_phi = phi
    R_z = np.array(
        [
            [np.cos(vinkel_phi), -np.sin(vinkel_phi), 0],
            [np.sin(vinkel_phi), np.cos(vinkel_phi), 0],
            [0, 0, 1]
        ], dtype=np.float64)

    # Rotation i y-led (theta).
    # För att z-axeln ska sammanfalla med riktningsvektorn
    # -> måste rotationsvinkeln vara theta_A.
    R_y = np.array(
        [
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)],
        ], dtype=np.float64)

    # Först rotation i y-led, sedan rotation i z-led.
    # R = np.dot(R_y, R_z)
    R = np.dot(R_z, R_y)

    return R
