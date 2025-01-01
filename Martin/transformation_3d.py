from imports import *


def transform(x_B, y_B, z_B, d_x, d_y, d_z, phi, theta):
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
    :param x_B: x-komponent av ny riktningsvektor (steg) i B's koordinatsystem
    :param y_B: y-komponent --''--
    :param z_B: z-komponent --''--
    :param d_x: x-komp av steg från A till B
    :param d_y: y-komp --''--
    :param d_z: z-komp --''--
    :param phi: vinkel sfärisk
    :param theta: vinkel sfärisk
    """
    enhets_vektorer_A = np.array([[1,0,0,0],[0,1,0,0], [0,0,1,0],[0,0,0,0]])

    # transformera först med rotation i z-led (phi)
    # 4D matrix för att kunna skapa homogenitet-transformsmatris
    R_z = np.array([[np.cos(phi), -np.sin(phi), 0, 0], [np.sin(phi), np.cos(phi), 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]])

    # transformera sedan i theta
    # tänk att theta rotationsmatrisen har anti-clockwise rotation som standard, men theta är clockwise, så ta -theta
    R_y = np.array(
        [[np.cos(theta), 0, -np.sin(theta), 0], [0, 1, 0, 0], [np.sin(theta), 0, np.cos(theta), 0], [0, 0, 0, 0]])

    # tänk att om theta > 0 måste kompenseringen vara rotation med -theta
    # R_x = np.array(
    #     [[1,0,0,0], [0, np.cos(theta), -np.sin(theta),0], [0, np.sin(theta), np.cos(theta), 0], [0, 0, 0, 0]])

    # R_y = np.identity(4)
    # R_z = np.identity(4)

    # R = R_y @ R_z
    R = R_z @ R_y
    # R = R_x @ R_z
    Homogenous_matrix = R + np.array([[0, 0, 0, d_x], [0, 0, 0, d_y], [0, 0, 0, d_z], [0, 0, 0, 1]])
    result = Homogenous_matrix @ np.array([[x_B], [y_B], [z_B], [1]])
    enhets_vektorer_B = Homogenous_matrix @ enhets_vektorer_A
    return result, enhets_vektorer_B


theta = pi / 8
phi = pi / 8
x_B = 1
y_B = 1
z_B = 1
s = 3
d_x = s * np.cos(theta) * np.cos(phi)
d_y = s * np.cos(theta) * np.sin(phi)
d_z = s * np.sin(theta)
print(np.cos(pi))




result, enhets_vektorer_B = transform(x_B, y_B, z_B, d_x, d_y, d_z, phi, theta)
x_A = result[0]
y_A = result[1]
z_A = result[2]

print(enhets_vektorer_B)
print(0.5**2 + 0.5**2 + 0.707**2)
print(0.707**2 + 0.707**2)

fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
ax.set_xlim([0, 4])
ax.set_ylim([0, 4])
ax.set_zlim([0, 4])

ax.quiver(0, 0, 0, 1, 0, 0, color='orange')
ax.quiver(0, 0, 0, 0, 1, 0, color='green')
ax.quiver(0, 0, 0, 0, 0, 1, color='green')

# ax.quiver(d_x, d_y, d_z, enhets_vektorer_B[0,0], enhets_vektorer_B[1,0], enhets_vektorer_B[2,0], color='green')
# ax.quiver(d_x, d_y, d_z, enhets_vektorer_B[0,1], enhets_vektorer_B[1,1], enhets_vektorer_B[2,1], color='green')
# ax.quiver(d_x, d_y, d_z, enhets_vektorer_B[0,2], enhets_vektorer_B[1,2], enhets_vektorer_B[2,2], color='green')

ax.quiver(d_x, d_y, d_z, enhets_vektorer_B[0,0], enhets_vektorer_B[1,0], enhets_vektorer_B[2,0], color='orange')
ax.quiver(d_x, d_y, d_z, enhets_vektorer_B[0,1], enhets_vektorer_B[1,1], enhets_vektorer_B[2,1], color='green')
ax.quiver(d_x, d_y, d_z, enhets_vektorer_B[0,2], enhets_vektorer_B[1,2], enhets_vektorer_B[2,2], color='green')

# ax.quiver(x_new, y_new, z_new, x_bar_new[0], x_bar_new[1], x_bar_new[2], color='red')
# ax.quiver(x_new, y_new, z_new, y_bar_new[0], y_bar_new[1], y_bar_new[2], color='red')
# ax.quiver(x_new, y_new, z_new, z_bar_new[0], z_bar_new[1], z_bar_new[2], color='red')

ax.quiver(0,0,0,x_A,y_A,z_A, color='red')
ax.quiver(0,0,0,d_x,d_y,d_z, color='blue')
ax.quiver(d_x, d_y, d_z, (x_A-d_x), (y_A-d_y), (z_A-d_z), color='magenta')


ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()

# x och z
print(np.dot(enhets_vektorer_B[0:2,0], enhets_vektorer_B[0:2,2]))
# x och y
print(np.dot(enhets_vektorer_B[0:2,0], enhets_vektorer_B[0:2,1]))
# y och z
print(np.dot(enhets_vektorer_B[0:2,1], enhets_vektorer_B[0:2,2]))