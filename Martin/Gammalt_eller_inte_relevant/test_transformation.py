from imports import *
from upg12_rotation_matris import rotations_matris
from upg12_steg_transformation import transformera_koordinatsystem
from upg12_förflyttning import förflyttning

iterationer = 50

x_list, y_list, z_list, dos_list = [], [], [], []

position_vektor = np.array([0, 0, 0])
x, y, z = position_vektor[0], position_vektor[1], position_vektor[2]

x_list.append(x)
y_list.append(y)
z_list.append(z)
dos_list.append(1)

theta = pi/25
phi = 0

steg = 1

R = rotations_matris(phi, theta)

i = 0
while i < iterationer:
    i += 1
    dos_list.append(i)
    dx, dy, dz, R = transformera_koordinatsystem(steg, phi, theta, steg, phi, theta, R)

    x, y, z, _, _, _ = förflyttning(x, y, z, dx, dy, dz)
    position_vektor = np.array([x, y, z])

    x_list.append(x)
    y_list.append(y)
    z_list.append(z)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
ax.scatter(x_list, y_list, z_list, c=dos_list, cmap='plasma', label='Partikel position')
# Fixa colorbar för att se energideponeringen i figuren

# fig.colorbar(ax=ax, label='Energideponering',)

ax.set_xlabel('x-axel (m)')
ax.set_ylabel('y-axel (m)')
ax.set_zlabel('z-axel (m)')

ax.legend()

# Visa figur
plt.tight_layout()
plt.show()
