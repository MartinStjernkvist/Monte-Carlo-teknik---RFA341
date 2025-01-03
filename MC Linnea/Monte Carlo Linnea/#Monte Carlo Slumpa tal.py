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






#Läser en Excel fil
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import openpyxl as op

#Hittar Excel filen i datorn och läs in
file_attenuering=pd.read_excel(r"C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Attenueringsdata.xlsx")
#Orelevanta datan som ska ta bort i filen som t.ex. titanium och bly

#Testar och se om det blir en förändring
print(file_attenuering.shape)

#Testa drop funktionen och ta bort 
file_attenuering=file_attenuering.drop([20, 21, 22, 24,25], axis=1) # Nope..
 
#Testar och se om det blir en förändring
print(file_attenuering.shape)


Foton_energi=208.3 #Test, ta bort sen,  energi på Lu 177 i keV
Energi_data=file_attenuering['E'] #Data i keV också


data_kropp=file_attenuering.values.tolist() #Resten av attenueringskoefficienten i fantomen på rader

#Tar ut attenueringskoefficienterna från Excel datan 
for i in range(len(Energi_data)):
    attenueringskoefficienter=[]
    if Energi_data[i]== int(Foton_energi):
        #Interpolera energi så man får samma värde från datan
        Energi_delar=np.linspace(Energi_data[i], Energi_data[i+1],5) 
        #Interpolera raden i och i+1 för alla datavärden
        data_kropp_delar=np.linspace(data_kropp[i],data_kropp[i+1],5)

        for j in range(len(Energi_delar)):
            if Energi_delar[j]==Foton_energi:
                #Sampla attenueringskoefficienten från alla vävnader i fantomkroppen 
                attenueringskoefficienter.append(data_kropp_delar[j]) #Hitta ett sätt att ta ta bort delar som inte är relevanta i filen
                

                break

    else:
        continue



import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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

olika_energier=np.linspace(np.min(Energi_Y90), np.max(Energi_Y90))

plt.plot(olika_energier,polynom_funktion(olika_energier,*params))

#plt.show()

#Inver transformera funktionen


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