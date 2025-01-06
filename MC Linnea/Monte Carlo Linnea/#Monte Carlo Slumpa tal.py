#Monte Carlo Slumpa tal
import random as rand
rand.random() # ger slumptal mellan 0 till 1
import numpy as np
np.random.rand( )# ger slumptal mellan 0 till 1)
import math as ma



import pandas as pd
#Testar att filerna fungerar:
file_tvärsnitt=pd.read_excel(r"C:\Users\Admin\Documents\GitHub\Monte-Carlo-teknik---RFA341\given_data\Tvärsnittstabell_Elektroner.xlsx")

Energi_data=file_tvärsnitt['Energy(eV)'].to_list()
#nope..


"""
#Artikel från Persliden, 1983 
#Vinkeln i Coherentspridning
while True:
    g_theata=np.random.random()#eller np.random random(0.5,1)?? 
    p=np.random.random()
    #Sampla x^2 från funktion f(x^2,Z) och får theata
    if g_theata>p:
        #Acceptera theata
        break
    else:
        continue
"""