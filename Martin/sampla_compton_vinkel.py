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
import math as ma


# foton_energin= 1 # Test energi på inkommande fotonen
def Compton_vinkel_och_energiförlust(foton_energi):
    while True:
        # Slumpa tre slumpmässiga tal
        R_1 = rand.random()
        R_2 = rand.random()
        R_3 = rand.random()

        # Khans metod (från artikeln "A Monte Carlo program for the simulation of scintillation camera characteristics") för att ta reda på vinkeln och energin på fotonen:
        if R_1 <= (2 * foton_energi + 1) / (2 * foton_energi + 9):
            
            nu = 2 * foton_energi * R_2
            if R_3 <= 4 * (nu ** -1 - nu ** -2):
                theta = ma.acos(1 - 2 * R_2)
                spridd_foton_energi = foton_energi / nu
                energi_förlust = foton_energi - spridd_foton_energi

                break

        elif R_1 > (2 * foton_energi + 1) / (2 * foton_energi + 9):

            nu = (2 * foton_energi + 1) / (2 * R_2 * foton_energi + 1)
            theta = ma.acos(1 - (nu - 1) / foton_energi)

            if R_3 <= 0.5 * (ma.cos(theta) ** 2 + nu ** -1):

                spridd_foton_energi = foton_energi / nu
                theta = ma.acos(1 - (nu - 1) / foton_energi)
                energi_förlust = foton_energi - spridd_foton_energi

                break
        else:
            continue

    return theta, spridd_foton_energi, energi_förlust
