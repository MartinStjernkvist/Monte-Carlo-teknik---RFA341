#Monte Carlo Slumpa tal

from imports import *

#from Elektron_energi import Närmast #VARFÖÖÖÖR???!!!! AHHHHHHHHHHHHHH


from Martin.upg1_sampla_energi_start import energi_start #Kolla varför det inte går att importera in
from Martin.upg1_matriser import slicad_fantom_matris, slicad_njure_matris, slicad_benmärg_matris
from Martin.upg1_sampla_position_start import position_start
from Martin.upg1_sampla_riktning_och_steg_start import riktning_uniform, steg
from Martin.upg1_sampla_steglängd import medelvägslängd
from Martin.upg1_sampla_växelverkan import växelverkan
from Martin.upg1_sampla_foto_vxv import foto_vxv
from Martin.upg1_bestäm_om_attenuerad import bestäm_om_attenuerad
from Martin.upg12_förflyttning import förflyttning
from Martin.upg1_attenueringsdata import attenueringsdata

df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)
df_anatomidefinitioner = pd.read_excel(anatomidefinitioner_file, index_col=None)
foton_energi = energi_start(Lu177_energi , Lu177_sannolikhet )
x_start, y_start, z_start = position_start(slicad_njure_matris)
voxel_värde = slicad_fantom_matris[x_start, y_start, z_start]
instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
mu = instans.mu_max()
theta, phi = riktning_uniform()
steglängd = medelvägslängd(mu)

dx, dy, dz = steg(theta, phi, steglängd)
x, y, z, x_round, y_round, z_round = förflyttning(x_start, y_start, z_start, dx, dy, dz)

voxel_värde=slicad_fantom_matris[x_round,y_round,z_round]
#Ta reda på vilken organ det är och attenueringskoefficienten my_voxel
mu_voxel=attenueringsdata.mu()

#Fråga om mu_max är i riktningen eller matrisen!, sätt sedan funktionen
mu_max=2
#Ta reda på om det sker en vävelverkan eller inte efter steglängden
while True:
    R=np.random.random()

    if R<= mu_voxel/mu_max:#mu_max är största voxelvärde i riktningen?
        #Växelverkan sker i steglängden och riktningen !!!
        print(steglängd) #return steglängd eller sätt funktion där man bestämmer vxv
        break
        
        
    else:
        continue #Sampla ny stegläng men i samma riktning





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