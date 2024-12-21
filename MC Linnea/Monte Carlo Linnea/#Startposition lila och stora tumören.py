#Startposition lila och stora tumören

import random as rand
rand.random() # ger slumptal mellan 0 till 1
import numpy as np
np.random.rand( )# ger slumptal mellan 0 till 1)
import math as ma
#Slumpa startpositionen för alfa och beta
#Antar att tumören är i mitten av koordinatsystemet

#Position på sfärens yta för lilla tumören
radie_liten=300*10**-6 # i meter
Sfäryta_liten=4*radie_liten**2*ma.pi
print(Sfäryta_liten)

Antal_iterationer=100
Startpos_yta=[]
for i in range(Antal_iterationer):
    x=np.random.uniform(-0.5,0.5)#Slumpmässig koordinat mellan 1 m i koordinatsystemet
    y=np.random.uniform(-0.5,0.5)
    z=np.random.uniform(-0.5,0.5)
    if radie_liten**2==y**2+z**2+x**2:
        Startpos_yta.append([x,y,z])
    else:
        continue
#Position i sfärens volym för stora tumören
radie_stor=1*10**-3 # i meter
Sfärvolym_stor=4*radie_liten**3*ma.pi/3
print(Sfärvolym_stor)

Antal_iterationer=100
Startpos_volym=[]
for i in range(Antal_iterationer):
    x=np.random.uniform(-0.5,0.5) #Slumpmässig koordinat mellan 1 m i koordinatsystemet
    y=np.random.uniform(-0.5,0.5)
    z=np.random.uniform(-0.5,0.5)
    if radie_stor**2>=y**2+z**2+x**2:
        Startpos_volym.append([x,y,z])
    else:
        continue
