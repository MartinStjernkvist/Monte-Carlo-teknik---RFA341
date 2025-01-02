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
Energi= file['Energy (MeV)'] #MeV
Intensitet=file['#/nt']
#Plottar ut värdena
plt.figure(1)
plt.scatter(Energi, Intensitet)
#Visa figuren

#Plottar ut punkterna i excelfilen och gör en kurvanpassning

def polynom_funktion(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d
params, cv= curve_fit (polynom_funktion,Energi, Intensitet)

olika_energier=np.linspace(np.min(Energi), np.max(Energi))

plt.plot(olika_energier,polynom_funktion(olika_energier,*params))

plt.show()

#Inver transformera funktionen