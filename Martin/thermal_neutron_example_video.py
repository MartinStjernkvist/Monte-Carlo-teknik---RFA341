import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import time
from scipy.interpolate import interp1d

'''video: https://www.youtube.com/watch?v=nVvFQ0e8J8g'''

start = time.time()

pi = np.pi

lambda_a = 45  # mean free path absorption
lambda_s = 0.3  # # mean free path scatter

sigma_a = 1 / lambda_a
sigma_s = 1 / lambda_s

sigma_t = sigma_a + sigma_s

lambda_t = 1 / sigma_t


n_particles = 10**3
x = np.zeros(n_particles)
y = np.zeros(n_particles)
z = np.zeros(n_particles)
r = np.zeros(n_particles)

for j in range(n_particles):

    # origin particle, should not be isotropic
    x[j] = 0
    y[j] = 0
    z[j] = 0

    # starting values
    is_absorbed = 0
    i = 0

    while is_absorbed == 0:

        s = - lambda_t * np.log(random.rand())  # sample

        # angles
        theta = np.arcsin(-1 + 2 * random.rand())
        phi = 2 * pi * random.rand()

        # movement
        dx = s * np.cos(theta) * np.cos(phi)
        dy = s * np.cos(theta) * np.sin(phi)
        dz = s * np.sin(theta)

        # new position
        x[j] = x[j] + dx
        y[j] = y[j] + dy
        z[j] = z[j] + dz

        i += 1

        # what happens at each step: absorbed or not?
        if random.rand() < sigma_a / sigma_t:
            is_absorbed = 1

    # distance traveled
    r[j] = np.sqrt(x[j]**2 + y[j]**2 + z[j]**2)

r_avg = np.mean(r)
print(r_avg)

'''
PLOTTING
'''

big_font_size = 25
small_font_size = big_font_size * 0.75

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
ax.scatter(x, y, z, color='red')
#ax.scatter(x, y, z, rstride=1, cstride=1, cmap='binary', edgecolor='none', label='no correction')
# ax.contour3D(V, AR, L_D_no_corr, 100, cmap='viridis')
# ax.plot_surface(V, AR, L_D_no_corr, rstride=1, cstride=1, cmap='viridis', edgecolor='none', label='corrected with RE')

ax.set_xlabel('x', fontsize=big_font_size)
ax.set_ylabel('y', fontsize=big_font_size)
ax.set_zlabel('z', fontsize=big_font_size)

# ax.set_zlim(0, 150)

plt.xticks(fontsize=small_font_size)
plt.yticks(fontsize=small_font_size)
plt.legend(fontsize=small_font_size)
plt.show()




end = time.time()

runtime = round((end - start), 1)
if runtime < 60:
    print(f'Runtime: {runtime} seconds')
else:
    print('Runtime: ' + str(round((runtime/60), 1)) + ' minutes')