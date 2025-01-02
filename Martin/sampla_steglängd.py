from imports import *


def invers_funktion(x, mu):
    # x borde i det här fallet vara rand(0,1)
    return -np.log(x) / mu


def medelvägslängd(mu):
    medelvägslängd = invers_funktion(np.random.rand(), mu)
    return medelvägslängd