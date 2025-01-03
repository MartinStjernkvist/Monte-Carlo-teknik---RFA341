#Monte Carlo Slumpa tal
import random as rand
rand.random() # ger slumptal mellan 0 till 1
import numpy as np
np.random.rand( )# ger slumptal mellan 0 till 1)
import math as ma



#Slumpa ut energi på betakällan
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import integrate
#Hittar filen och ta info från excelfilen
#Läser en Excel fil
file_Y90=pd.read_excel(r"C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Y90_Spektrum.xlsx" )
#För att se att tabellen stämmer:
# print(file)

#Ta ut värderna på energin och intensiteten betakällan
Energi_Y90= file_Y90['Energy (MeV)'] #MeV
Intensitet_Y90=file_Y90['#/nt']
#Plottar ut värdena
plt.figure(1)
plt.scatter(Energi_Y90, Intensitet_Y90)
#Visa figuren

#Anta att tumören är i en vävnad alltså se på vatten för tvärsnitten

#Plottar ut punkterna i excelfilen och gör en kurvanpassning

def polynom_funktion(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d
params, cv= curve_fit (polynom_funktion,Energi_Y90, Intensitet_Y90)
print(*params)
olika_energier=np.linspace(np.min(Energi_Y90), np.max(Energi_Y90))

plt.plot(olika_energier,polynom_funktion(olika_energier,*params))

plt.show()

Skärpunkt=2.206882599192 #Enligt Wolfram alfa

def Lösning_tredjegradare(a,b,c,d):
    p=c/a-((b/a)**2)/3
    q=d/a+(2*(b/a)**3-9*b*c/a**2)/27
    D=(p/3)**3+(q/2)**2
    u=(-q/2+D**0.5)**(1/3)
    v=(-q/2-D**0.5)**(1/3)
    x=u+v-b/(3*a)
    return x
print(Lösning_tredjegradare(*params)) #Stämmer inte med grafen

#Inver transformera funktionen
#k=1/integrate.quad(polynom_funktion(olika_energier,*params),0,Skärpunkt)
#print(k)

"""
#Sfärisk koordinat för riktningnen som fotoner kan färdas i
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
steg=1 #Test steglängd
Antal_iterationer=1000

fig = plt.figure(2)
ax = plt.axes(projection='3d')

#Testa från samplingskoden
#from Martin.sampla_riktning_och_steg_start import riktning_start

def riktning_start():
    # theta = random.gauss(pi / 2, 1)  # theta = np.arcsin(-1 + 2 * random.rand())

    theta = np.arccos(-1 + 2 * rand.random()) 
    phi = 2 * ma.pi * rand.random()
    return theta, phi

for i in range(Antal_iterationer):     
    #Slumpar vinklarna på teata och si
    theata,phi=riktning_start()
    #theata=rand.gauss(ma.pi/2,1) #Fortfarande många vid toppen och botten när många iterationer görs, men är bättre än uniform
    #phi=np.random.uniform(0,2*ma.pi) #Slumpar uniform phi
    #Koordinaterna för vinklarna

    x=steg*ma.sin(theata)*ma.cos(phi)
    y=steg*ma.sin(theata)*ma.sin(phi)
    z=steg*ma.cos(theata)
    #Testar för att se hur punkerna är fördelade
    ax.scatter(x,y,z,color='blue',s=2)#Plottar ut punkter
   

#Plotta figur
plt.show()
"""