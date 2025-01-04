from imports import *

# #Sampla Compton vinkel och den sprida fotonenergin
# import random as rand
# rand.random() # ger slumptal mellan 0 till 1
# import numpy as np
# np.random.rand( )# ger slumptal mellan 0 till 1)
# import math as ma
#
# #foton_energin= 1 # Test energi på inkommande fotonen
# def Compton_vinkel_och_energiförlust(foton_energin):
#     while True:
#         #Slumpa tre slumpmässiga tal
#         R_1=rand.random()
#         R_2=rand.random()
#         R_3=rand.random()
#         #Khans metod (från artikeln "A Monte Carlo program for the simulation of scintillation camera characteristics") för att ta reda på vinkeln och energin på fotonen:
#         if R_1<=(2*foton_energin+1)/(2*foton_energin+9):
#             nu=2*foton_energin*R_2
#             if R_3<=4*(nu**-1-nu**-2):
#                 theata=ma.acos(1-2*R_2)
#                 sprid_foton_energin=foton_energin/nu
#                 energi_förlust=foton_energin-sprid_foton_energin
#                 break
#             return theata, sprid_foton_energin, energi_förlust
#
#         elif R_1> (2*foton_energin+1)/(2*foton_energin+9):
#             nu=(2*foton_energin+1)/(2*R_2*foton_energin+1)
#             theata=ma.acos(1-(nu-1)/foton_energin)
#
#             if R_3<=0.5*(ma.cos(theata)**2+nu**-1):
#                 sprid_foton_energin=foton_energin/nu
#                 theata=ma.acos(1-(nu-1)/foton_energin)
#                 energi_förlust=foton_energin-sprid_foton_energin
#                 break
#             return theata, sprid_foton_energin, energi_förlust
#         else:
#             continue


# Sampla Compton vinkel och den sprida fotonenergin
import random as rand


# foton_energin= 1 # Test energi på inkommande fotonen
@jit(nopython=True)
def compton_vinkel_och_energiförlust(foton_energi):
    """
    Khans rejektionsalgoritm, artikeln "A Monte Carlo program for the simulation of scintillation camera characteristics"
    Algoritmen erhåller vinkeln och energin på fotonen.
    """

    while True:
        # Slumpa tre slumpmässiga tal
        R_1 = np.random.rand()
        R_2 = np.random.rand()
        R_3 = np.random.rand()

        alpha = foton_energi / E_e


        if R_1 <= (2 * alpha + 1) / (2 * alpha + 9):

            eta = 2 * alpha * R_2
            if R_3 <= 4 * (eta ** -1 - eta ** -2):
                theta = np.arccos(1 - 2 * R_2)
                spridd_foton_energi = foton_energi / eta
                energi_förlust = foton_energi - spridd_foton_energi

                break

        elif R_1 > (2 * alpha + 1) / (2 * alpha + 9):

            eta = (2 * alpha + 1) / (2 * R_2 * alpha + 1)
            theta = np.arccos(1 - (eta - 1) / alpha)

            if R_3 <= 0.5 * (np.cos(theta) ** 2 + eta ** -1):

                spridd_foton_energi = foton_energi / eta
                theta = np.arccos(1 - (eta - 1) / alpha)
                energi_förlust = foton_energi - spridd_foton_energi

                break

        else:
            continue

    return theta, spridd_foton_energi, energi_förlust
