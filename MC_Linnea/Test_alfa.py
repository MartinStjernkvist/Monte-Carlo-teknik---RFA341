#Monte Carlo Slumpa tal

from imports import *
#Riktning är uniform för alfa, och kommer gå i samm riktning utan att ändra på vinkeln efter varje kollision
from Martin.upg1_sampla_riktning_och_steg_start import riktning_uniform # Går inte att få upp VARFÖÖÖÖR???!!!! AHHHHHHHHHHHHHH
from Martin.upg1_sampla_energi_start import energi_start

S=1*10**(-6) #1 um, hitta steglängden för energin 
steglängd=np.random.random()*S

Energi=energi_start(At211_energi,At211_sannolikhet) #Samplar på samma sätt alfa energin som fotonerna



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