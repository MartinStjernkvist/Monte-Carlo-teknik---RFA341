#Slumpmässiga Alfaenergier

#Kan annars anta att energin är bara 5869 keV p.g.a högst sannolikhet

import random as rand
rand.random() # ger slumptal mellan 0 till 1
import numpy as np
np.random.rand( )# ger slumptal mellan 0 till 1)
import math as ma

#Ursprungsenergier från nukliden
At211_energ=[5869, 5211.9,5140.3, 4993.4,4895.4] #keV

#Intensitet i % för energierna
At211_intens=[41.78, 0.0039, 0.0011, 0.0004, 0.00004]
#Gör om intensiteten till tal istället för %

for i in range(len(At211_intens)):
    At211_intens[i]=At211_intens[i]*10**-2

#Skapar en lista med sannolikheten att en viss energi på fotonen skapas
At211_sannolik=np.zeros(len(At211_intens))


for i in range(len(At211_intens)):
    if i==0:
        At211_sannolik[i]=At211_intens[i]
    else:
        At211_sannolik[i]= At211_intens[i]+At211_intens[i-1]

#Slumpa ut energin på alfapartikeln
Elektron_energi=[]
Antal_iterationer=100
for i in range(Antal_iterationer):
    Slump_tal=np.random.rand() #Slumptal
    if Slump_tal<=At211_sannolik[0]:
        Elektron_energi.append(At211_energ[0]) 
    elif Slump_tal<=At211_sannolik[1]:
        Elektron_energi.append(At211_energ[1]) 
    elif Slump_tal<=At211_sannolik[2]:
        Elektron_energi.append(At211_energ[2])
    elif Slump_tal<=At211_sannolik[3]:
        Elektron_energi.append(At211_energ[3]) 
    elif Slump_tal<=At211_sannolik[4]:
        Elektron_energi.append(At211_energ[4]) 
    else:
        continue
#print(Elektron_energi)