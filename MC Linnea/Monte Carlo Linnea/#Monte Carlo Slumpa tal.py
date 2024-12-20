#Monte Carlo Slumpa tal
import random as rand
rand.random() # ger slumptal mellan 0 till 1
import numpy as np
np.random.rand( )# ger slumptal mellan 0 till 1)
import math as ma

#Ursprungsenergier från nukliden
Y90_energ=2278.7 #keV men ska slumpa i spektrumet i en av filerna



#Steglängden på fotonerna

#Ta reda på majortiteta attenuerigskoefficienten för varje matris/slice
file_path_attenuering="C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Attenueringsdata.xlsx"
... #ta reda på my_m genom att ta max för allavoxlar attenueringskoefficient
my_m=np.max()

#Steglängden på fotonerna
steg=-ma.log(rand.random())/my_m


#Hittar filen 
#file_path_betasönderfall="C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Y90_Spektrum.xlsx"

#Läser en Excel fil
import pandas as pd
df=pd.read_excel(r"C:\Users\Admin\Documents\GitHub\Monte Carlo Linnea\Y90_Spektrum.xlsx" )
print(df)
#Lägg till sannolikheten för beta sönderfall? på något sätt hitta hur man gör det i en excel fil

#Position på sfärens yta för lilla tumören
radie_liten=300*10**-6
Sfär_liten=4*radie_liten**2*ma.pi

#position i sfärens volym för stora tumören
radie_stor=200*10**-2
Volym_stor=4*radie_liten**3*ma.pi/3
