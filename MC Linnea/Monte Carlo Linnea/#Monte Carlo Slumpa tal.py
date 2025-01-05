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
from scipy.special import cbrt 

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
#print(*params)
olika_energier=np.linspace(np.min(Energi_Y90), np.max(Energi_Y90))

plt.plot(olika_energier,polynom_funktion(olika_energier,*params))

#plt.show()

Skärpunkt=2.206882599192 #Enligt Wolfram alfa

def Lösning_tredjegradare(a,b,c,d):
    p=c/a-((b/a)**2)/3
    q=d/a+(2*(b/a)**3-9*b/a*c/a)/27
    D=(p/3)**3+(q/2)**2
    u=cbrt(-q/2+D**0.5)
    v=cbrt(-q/2-D**0.5)
    x=u+v-b/(3*a)
    return x
#print(Lösning_tredjegradare(*params)) #Stämmer inte med grafen

#Inver transformera funktionen
#k=1/integrate.quad(polynom_funktion(olika_energier,*params),0,Skärpunkt) #kunde inte integrera av någon anledning
#print(k)

#Testar att filerna fungerar:
file_tvärsnitt=pd.read_excel(r"C:\Users\Admin\Documents\GitHub\Monte-Carlo-teknik---RFA341\given_data\Tvärsnittstabell_Elektroner.xlsx")
Energi_data=file_tvärsnitt['Energy(eV)'].to_list()
    


