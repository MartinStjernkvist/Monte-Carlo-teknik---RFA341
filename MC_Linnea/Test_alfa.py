#Monte Carlo Slumpa tal

from imports import *
#Riktning är uniform för alfa, och kommer gå i samm riktning utan att ändra på vinkeln efter varje kollision
#from Martin.upg1_sampla_riktning_och_steg_start import riktning_uniform # Går inte att få upp VARFÖÖÖÖR???!!!! AHHHHHHHHHHHHHH
#from Martin.upg1_sampla_energi_start import energi_start

S=1*10**(-6) #1 um, hitta steglängden för energin i elektronen
steglängd=np.random.random()*S

#Energi=energi_start(At211_energi,At211_sannolikhet) #Samplar på samma sätt alfa energin som fotonerna, är i MeV

def position_start_alpha_innanför(radie_sfär, phi, theta):
    r = radie_sfär* np.random.rand()

    x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return x,y,z # position_vektor

def position_start_alpha_skal(radie_sfär, phi, theta):
    radie_alpha=1.2*10**(-15)*4**(1/3) #radie i meter
    r = radie_sfär - 0.5 * radie_alpha  # För att inte endast theta = pi ska ge utslag

    x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)
    y = 0
    z = 0

    #position_vektor = np.array([x, y, z])
    return x,y,z #position_vektor

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
fig = plt.figure(1)
ax = plt.axes(projection='3d')
Antal_iterationer=100

radie_stor=1*10**-3 # i meter
radie_liten=300*10**-6 # i meter

def riktning_alpha():
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi

theata,phi=riktning_alpha()

for i in range(Antal_iterationer):
    x,y,z=position_start_alpha_innanför(radie_stor,phi, theata)
    ax.scatter(x,y,z,color='blue', s=3)
fig2 = plt.figure(2)
ax2 = plt.axes(projection='3d') 
for i in range(Antal_iterationer):
    ax2.scatter(position_start_alpha_skal(radie_liten,phi,theata),color='green', s=3)


#Visa figur
plt.show()
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