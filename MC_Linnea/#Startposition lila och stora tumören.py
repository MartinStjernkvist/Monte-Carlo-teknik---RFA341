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


Antal_iterationer=10000
Startpos_yta=[]
#Vill plotta ut punkterna för att se att de stämmer med sfären
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
fig = plt.figure()
ax = plt.axes(projection='3d')

for i in range(Antal_iterationer):
    x=np.random.uniform(-0.5*10**-3,0.5*10**-3)#Slumpmässig koordinat mellan 1 mm i koordinatsystemet, värden i meter
    y=np.random.uniform(-0.5*10**-3,0.5*10**-3)
    z=np.random.uniform(-0.5*10**-3,0.5*10**-3)
    if radie_liten**2>=y**2+z**2+x**2 and y**2+z**2+x**2> (299*10**-6)**2 :
        Startpos_yta.append([x,y,z])
        #ax.scatter(x,y,z,color='green', s=3)
        continue
    else:
        continue
#print(Startpos_yta)


#Position i sfärens volym för stora tumören
radie_stor=1*10**-3 # i meter
Sfärvolym_stor=4*radie_stor**3*ma.pi/3

Antal_iterationer=100
Startpos_volym=[]
for i in range(Antal_iterationer):
    x=np.random.uniform(-1*10**-3,1*10**-3) #Slumpmässig koordinat mellan 2 mm i koordinatsystemet
    y=np.random.uniform(-1*10**-3,1*10**-3)
    z=np.random.uniform(-1*10**-3,1*10**-3)
    if radie_stor**2>=y**2+z**2+x**2:
        Startpos_volym.append([x,y,z])
        ax.scatter(x,y,z,color='blue', s=3) 
        continue
    else:
        continue
#print(Startpos_volym)
#Visa figur
plt.show()

def Startposition_liten_tumör():
    radie_liten=300*10**-6 # i meter
    while True:
        #Slumpar en position från -0.5 mm till 0.5 mm
        x=np.random.uniform(-0.5*10**-3,0.5*10**-3)#Slumpmässig koordinat mellan 1 mm i koordinatsystemet, värden i meter
        y=np.random.uniform(-0.5*10**-3,0.5*10**-3)
        z=np.random.uniform(-0.5*10**-3,0.5*10**-3)
        if radie_liten**2>=y**2+z**2+x**2 and y**2+z**2+x**2> (299*10**-6)**2 :
            break
    return x,y,z

def Startposition_stor_tumör():
    radie_stor=1*10**-3 # i meter
    while True:
        #Slumpar tal mellan -1 mm och 1 mm 
        x=np.random.uniform(-1*10**-3,1*10**-3) #Slumpmässig koordinat mellan 2 mm i koordinatsystemet
        y=np.random.uniform(-1*10**-3,1*10**-3)
        z=np.random.uniform(-1*10**-3,1*10**-3)
        if radie_stor**2>=y**2+z**2+x**2:
            break
    return x,y,z