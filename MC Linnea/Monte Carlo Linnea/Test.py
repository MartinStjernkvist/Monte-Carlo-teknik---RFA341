#Monte Carlo Slumpa tal

from imports import *
#from Elektron_energi import Närmast #VARFÖÖÖÖR???!!!! AHHHHHHHHHHHHHH

#Testar att filerna fungerar:
file_scatter=pd.read_excel(r'MC Linnea/Monte Carlo Linnea/Monte Carlo, test.xlsx',sheet_name='Scattering power') 


Energi_data=file_scatter['E(eV)'].to_list() #Hittar inte kolumnerna i xlsx filen
Scatterpower_list=file_scatter['Scatter power Water'].to_list() #T/rho för vatten


#kanske inte relevant att ha Rayleight spridning för det verkar som att majoriteten av vxv är antingen compton eller fotovxv
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