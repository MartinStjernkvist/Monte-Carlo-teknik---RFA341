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
df=pd.read_excel(r"C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Y90_Spektrum.xlsx" )
print(df)
#Plottar ut punkterna i excelfilen och gör en kurvanpassning


#Sfärisk koordinat för riktningnen som fotoner kan färdas i
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
steg=1 #Test steglängd
Antal_iterationer=500
Punkter=[]
fig = plt.figure()
ax = plt.axes(projection='3d')

for i in range(Antal_iterationer):     
    #Slumpar vinklarna på teata och si
    #theata=np.random.uniform(0,ma.pi) #fixa detta 
    theata=rand.gauss(ma.pi/2,1) #Fortfarande många vid toppen och botten när många iterationer görs
    si=np.random.uniform(0,2*ma.pi)
    #Koordinaterna för vinklarna
    x=steg*ma.sin(theata)*ma.cos(si)
    y=steg*ma.sin(theata)*ma.sin(si)
    z=steg*ma.cos(theata)
    ax.scatter(x,y,z,color='blue',s=2)#Plottar ut punkter
    #Punkter.append([x,y,z])


#Plotta figur    
plt.show()