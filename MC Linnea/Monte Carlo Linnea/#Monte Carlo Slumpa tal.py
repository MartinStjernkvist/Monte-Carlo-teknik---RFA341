#Monte Carlo Slumpa tal
import random as rand
rand.random() # ger slumptal mellan 0 till 1
import numpy as np
np.random.rand( )# ger slumptal mellan 0 till 1)
import math as ma

#Ursprungsenergier från nukliden
Y90_energ=2278.7 #keV men ska slumpa i spektrumet i en av filerna


#Steglängden på fotonerna, Wendy ska göra detta

#Ta reda på majortiteta attenuerigskoefficienten för varje matris/slice
#file_path_attenuering="C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Attenueringsdata.xlsx"
... #ta reda på my_m genom att ta max för allavoxlar attenueringskoefficient
#my_m=np.max()

#Steglängden på fotonerna
#steg=-ma.log(rand.random())/my_m




#Hittar filen och ta info från excelfilen
#file_path_betasönderfall="C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Y90_Spektrum.xlsx"

#Läser en Excel fil
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#Hittar Excel filen i datorn
file=pd.read_excel(r"C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Y90_Spektrum.xlsx" )
#För att se att tabellen stämmer:
# print(file)

#Ta ut värderna på energin och intensiteten
Energi= file['Energy (MeV)']
Intensitet=file['#/nt']
#Plottar ut värdena
plt.figure(1)
plt.scatter(Energi, Intensitet)
#Visa figuren
plt.show()
#Plottar ut punkterna i excelfilen och gör en kurvanpassning

def polynom_funktion(x,a,b,c):
    return a*x**2+b*x+c
params, cv= curve_fit (polynom_funktion,Energi, Intensitet)

olika_energier=np.linspace(np.min(Energi), np.max(Energi))


#Sampla compton vinkel och energiförlusten med Kahns metod i artikeln "A Monte Carlo program for the simulation of scintillation camera characteristics."

foton_energin= 1 # Test energi på inkommande fotonen
def Compton_vinkel_och_energiförlust():
    while True:
        #Slumpa tre slumpmässiga tal
        R_1=rand.random()
        R_2=rand.random()
        R_3=rand.random()
        if R_1<=(2*foton_energin+1)/(2*foton_energin+9):
            nu=2*foton_energin*R_2
            if R_3<=4*(nu**-1-nu**-2):
                theata=ma.acos(1-2*R_2)
                sprid_foton_energin=foton_energin/nu
                energi_förlust=foton_energin-sprid_foton_energin
                break
        else:
            continue
        if R_1> (2*foton_energin+1)/(2*foton_energin+9):
            nu=(2*foton_energin+1)/(2*R_2*foton_energin+1)
            theata=ma.acos(1-(nu-1)/foton_energin)
            if R_3<=0.5*(ma.cos(theata)**2+nu**-1):
                sprid_foton_energin=foton_energin/nu
                theata=ma.acos(1-(nu-1)/foton_energin)
                energi_förlust=foton_energin-sprid_foton_energin
                break
        else:
            continue




"""
#Sfärisk koordinat för riktningnen som fotoner kan färdas i
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
steg=1 #Test steglängd
Antal_iterationer=500

fig = plt.figure()
ax = plt.axes(projection='3d')


#Testa från samplingskoden

def riktning_start(steglängd): 
    theta = np.arcsin(-1 + 2 * rand.random())# Hitta ett sätt att sampla i alla riktningar för theata och lika fördelad, blir en parabol
    phi = 2 * ma.pi * rand.random()
    return theta, phi


for i in range(Antal_iterationer):     
    #Slumpar vinklarna på teata och si
    
    theata=rand.gauss(ma.pi/2,1) #Fortfarande många vid toppen och botten när många iterationer görs, men är bättre än uniform
    phi=np.random.uniform(0,2*ma.pi) #Slumpar uniform phi
    #Koordinaterna för vinklarna

    x=steg*ma.sin(theata)*ma.cos(phi)
    y=steg*ma.sin(theata)*ma.sin(phi)
    z=steg*ma.cos(theata)
    #Testar för att se hur punkerna är fördelade
    ax.scatter(x,y,z,color='blue',s=2)#Plottar ut punkter
   

#Plotta figur    
plt.show()
"""