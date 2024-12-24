from imports import *

def steg(theta, phi, steglängd, x,y,z):
    dx = steglängd * np.cos(theta) * np.cos(phi)
    dy = steglängd * np.cos(theta) * np.sin(phi)
    dz = steglängd * np.sin(theta)

    # new position
    x = x + dx
    y = y + dy
    z = z + dz
    return x, y, z