from imports import *
from upg12_rotation_matris import rotations_matris


def plotta_transformation(steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B, R_A):
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
    # H[:3, 3] = np.array([x_B_A, y_B_A, z_B_A], dtype=np.float64)
    H[:3, 3] = np.array([vektor_B_A[0], vektor_B_A[1], vektor_B_A[2]], dtype=np.float64)

    #   -----------------------------------
    #   Vektorer.
    #   -----------------------------------

    # vektor_C (A)
    vektor_C_A = np.dot(H, np.array([vektor_C_B[0], vektor_C_B[1], vektor_C_B[2], 1.0], dtype=np.float64))

    # vektor_B (A)
    # vektor_B_A = np.array([x_B_A, y_B_A, z_B_A, 1.0], dtype=np.float64)

    # vektor_B_A fast med en etta i slutet.
    vektor_B_A = np.array([vektor_B_A[0], vektor_B_A[1], vektor_B_A[2], 1.0])

    # Vill ha vektor_BC (A).
    # Vektorn från B till C, i A's koordinatsystem.
    vektor_BC_A = vektor_C_A - vektor_B_A

    # Vektorkomponenter av vektor_BC (A).
    x_BC_A, y_BC_A, z_BC_A = vektor_BC_A[0], vektor_BC_A[1], vektor_BC_A[2]

    # Genom att dela upp i vektorkomponenter kan positionen
    # för partikeln följas i det övegripande koordinatsystemet,
    # genom att addera stegvektorer.
    dx, dy, dz = x_BC_A, y_BC_A, z_BC_A

    return dx, dy, dz, R, vektor_B_A_innan, vektor_B_A, vektor_C_A, vektor_C_B, vektor_BC_A, R_innan, R_A


#   -----------------------------------
#
#   -----------------------------------

föregående_vektor = np.array([1, 2, 1])
vektor_A = np.array([2, 1, 2])

steg_A_B = 2
phi_A = 5 * pi / 8
theta_A = pi / 4

steg_B_C = 3
phi_B = pi / 4
theta_B = 2 * pi / 4

R_A = rotations_matris(-pi / 4, pi / 4)

eye = np.eye(3)
eye_R_A = np.dot(R_A, eye)

dx, dy, dz, R, vektor_B_A_innan, vektor_B_A, vektor_C_A, vektor_C_B, vektor_BC_A, R_innan, R_A = plotta_transformation(
    steg_A_B, phi_A, theta_A, steg_B_C, phi_B, theta_B, R_A)

eye_R_B = np.dot(R, eye)

font_size = 20

x_lim = 4
y_lim = x_lim
z_lim = x_lim


azim, elev = -115, 30
# -40, 30


#   -----------------------------------
#   Fig: A + rotation i A.
#   -----------------------------------


fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
ax.set_xlim([0, x_lim])
ax.set_ylim([0, y_lim])
ax.set_zlim([0, z_lim])

ax.scatter(föregående_vektor[0], föregående_vektor[1], föregående_vektor[2], label='Föregående punkt', color='deeppink',
           s=100)
ax.scatter(vektor_A[0], vektor_A[1], vektor_A[2], label='A', color='blue', s=100)
ax.scatter(0, 0, 0, label='Origo', color='black', s=100)

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 1, 0, 0, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 0, 1, 0, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 0, 0, 1, color='green', label='z-axel, A\'s koord-syst')

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 0], eye_R_A[1, 0], eye_R_A[2, 0], alpha=0.35, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 1], eye_R_A[1, 1], eye_R_A[2, 1], alpha=0.35, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 2], eye_R_A[1, 2], eye_R_A[2, 2], alpha=0.35, color='green',
          label='z-axel, A\'s virtuella koord-syst')

ax.quiver(0, 0, 0, 1, 0, 0, color='black')
ax.quiver(0, 0, 0, 0, 1, 0, color='black')
ax.quiver(0, 0, 0, 0, 0, 1, color='black')

ax.quiver(föregående_vektor[0], föregående_vektor[1], föregående_vektor[2], vektor_A[0] - föregående_vektor[0],
          vektor_A[1] - föregående_vektor[1], vektor_A[2] - föregående_vektor[2], color='darkviolet',
          label='Föregående steg -> A')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.azim, ax.elev = azim, elev

plt.title('Punkt A', fontsize=font_size * 1.5)
plt.tight_layout()
ax.legend(fontsize=font_size)
plt.savefig('Punkt_A')
plt.show()

#   -----------------------------------
#   Fig: Steg 1.
#   -----------------------------------

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
ax.set_xlim([0, x_lim])
ax.set_ylim([0, y_lim])
ax.set_zlim([0, z_lim])

ax.scatter(vektor_A[0], vektor_A[1], vektor_A[2], label='A', color='blue', s=100)
ax.scatter(0, 0, 0, label='Origo', color='black', s=100)

ax.scatter(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2], label='B',
           color='darkcyan', s=100)

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 1, 0, 0, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 0, 1, 0, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 0, 0, 1, color='green', label='z-axel, A\'s koord-syst')

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 0], eye_R_A[1, 0], eye_R_A[2, 0], alpha=0.35, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 1], eye_R_A[1, 1], eye_R_A[2, 1], alpha=0.35, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 2], eye_R_A[1, 2], eye_R_A[2, 2], alpha=0.35, color='green',
          label='z-axel, A\'s virtuella koord-syst')

ax.quiver(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2], eye_R_B[0, 0],
          eye_R_B[1, 0], eye_R_B[2, 0], color='sienna')
ax.quiver(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2], eye_R_B[0, 1],
          eye_R_B[1, 1], eye_R_B[2, 1], color='sienna')
ax.quiver(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2], eye_R_B[0, 2],
          eye_R_B[1, 2], eye_R_B[2, 2], color='orange',
          label='z-axel, B\'s koord-syst')

ax.quiver(0, 0, 0, 1, 0, 0, color='black')
ax.quiver(0, 0, 0, 0, 1, 0, color='black')
ax.quiver(0, 0, 0, 0, 0, 1, color='black')

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], vektor_B_A[0], vektor_B_A[1], vektor_B_A[2], color='olive',
          label='steg A -> B i A\'s koord-syst')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.azim, ax.elev = azim, elev

plt.title('Steg A -> B', fontsize=font_size * 1.5)
plt.tight_layout()
ax.legend(fontsize=font_size)
plt.savefig('Steg_A_B')
plt.show()

#   -----------------------------------
#   Fig: Steg 2.
#   -----------------------------------

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
ax.set_xlim([0, x_lim])
ax.set_ylim([0, y_lim])
ax.set_zlim([0, z_lim])

ax.scatter(vektor_A[0], vektor_A[1], vektor_A[2], label='A', color='blue', s=100)
ax.scatter(0, 0, 0, label='Origo', color='black', s=100)

ax.scatter(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2], label='B',
           color='darkcyan', s=100)
ax.scatter(vektor_A[0] + vektor_C_A[0], vektor_A[1] + vektor_C_A[1], vektor_A[2] + vektor_C_A[2], label='C',
           color='slategrey', s=100)

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 1, 0, 0, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 0, 1, 0, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], 0, 0, 1, color='green', label='z-axel, A\'s koord-syst')

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 0], eye_R_A[1, 0], eye_R_A[2, 0], alpha=0.35, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 1], eye_R_A[1, 1], eye_R_A[2, 1], alpha=0.35, color='red')
ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], eye_R_A[0, 2], eye_R_A[1, 2], eye_R_A[2, 2], alpha=0.35, color='green',
          label='z-axel, A\'s virtuella koord-syst')

ax.quiver(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2], eye_R_B[0, 0],
          eye_R_B[1, 0], eye_R_B[2, 0], color='sienna')
ax.quiver(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2], eye_R_B[0, 1],
          eye_R_B[1, 1], eye_R_B[2, 1], color='sienna')
ax.quiver(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2], eye_R_B[0, 2],
          eye_R_B[1, 2], eye_R_B[2, 2], color='orange',
          label='z-axel, B\'s koord-syst')

ax.quiver(0, 0, 0, 1, 0, 0, color='black')
ax.quiver(0, 0, 0, 0, 1, 0, color='black')
ax.quiver(0, 0, 0, 0, 0, 1, color='black')

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], vektor_B_A[0], vektor_B_A[1], vektor_B_A[2], color='olive',
          label='steg A -> B i A\'s koord-syst')

ax.quiver(vektor_A[0], vektor_A[1], vektor_A[2], vektor_C_A[0], vektor_C_A[1], vektor_C_A[2], color='magenta',
          label='steg A -> C i A\'s koord-syst')

ax.quiver(vektor_A[0] + vektor_B_A[0], vektor_A[1] + vektor_B_A[1], vektor_A[2] + vektor_B_A[2],
          vektor_C_A[0] - vektor_B_A[0], vektor_C_A[1] - vektor_B_A[1], vektor_C_A[2] - vektor_B_A[2],
          color='mediumpurple',
          label='steg B -> C i A\'s koord-syst')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.azim, ax.elev = azim, elev

plt.title('Steg B -> C', fontsize=font_size * 1.5)
plt.tight_layout()
ax.legend(fontsize=font_size)
plt.savefig('Steg_B_C')
plt.show()
