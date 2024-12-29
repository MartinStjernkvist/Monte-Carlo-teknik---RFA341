from imports import *

def riktning_start(steglängd):
    theta=random.gauss(pi/2,1) #theta = np.arcsin(-1 + 2 * random.rand())
    phi = 2 * pi * random.rand()
    return theta, phi

