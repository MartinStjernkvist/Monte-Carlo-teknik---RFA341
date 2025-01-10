from imports import *
#Riktning är uniform för alfa, och kommer gå i samm riktning utan att ändra på vinkeln efter varje kollision
#from Martin.upg1_sampla_riktning_och_steg_start import riktning_uniform # Går inte att få upp VARFÖÖÖÖR???!!!! AHHHHHHHHHHHHHH
#from Martin.upg1_sampla_energi_start import energi_start

#Energi=energi_start(At211_energi,At211_sannolikhet) #Samplar på samma sätt alfa energin som fotonerna, är i MeV
from Alfa_stp_och_RCSDA import Stopping_power_och_steglängd
from upg1_sampla_energi_start import energi_start

#Visa funktioner går att sätta jit men inte alla
@jit(nopython=True)
def riktning_alpha():
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi


#@jit(nopython=True)
def energiförlust_alpha(energi, steg):
    # Implementera stopping power
    # print('WIP')
    STP,_=Stopping_power_och_steglängd(energi)
    energiförlust=STP*steg #i MeV
    #energiförlust = energi * 0.1
    energi -= energiförlust

    if energi <= 0:
        energi = 0

    return energi


@jit(nopython=True)
def position_start_alpha_innanför(radie_sfär, phi, theta):
    r = radie_sfär * np.random.rand()

    x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
def position_start_alpha_skal(radie_sfär, phi, theta):
    r = radie_sfär - 0.5 * radie_alpha  # För att inte endast theta = pi ska ge utslag

    x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


#@jit(nopython=True)
def laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie, max_antal_steg=100):
    position_vektor = start_position
    energi = start_energi

    # trajectory = [tuple(position_vektor)]

    steg_storlek = steglängd / max_antal_steg

    riktning = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.cos(phi), np.cos(theta)])
    norm=np.sqrt(np.sum(riktning**2))
    riktning /= norm
    steg_vektor = riktning * steg_storlek


    # Under tiden som partikeln fortfarnade inte tagit hela sitt steg.
    for i in range(max_antal_steg):

        # print('steg_vektor', steg_vektor)
        position_vektor += steg_vektor
        #energi_förlust = energiförlust_alpha(energi, steg_storlek)
        energi = energiförlust_alpha(energi, steg_storlek)
       
       #Plottvärden
        energideponering_list=[] #för färgen av energideponering
        x=[]
        y=[]
        z=[]
        if np.dot(position_vektor, position_vektor) <= radie:
            innanför = True
            # trajectory.append(tuple(position_vektor))
            #print(f'Energideponering i position ', position_vektor)
            
            #Vill visa positionen för energideponeringen med en färg som motsvarar energideponeringen
            
            energideponering_list.append(energi-energiförlust_alpha(energi, steg_storlek))
            x.append(position_vektor[0])
            y.append(position_vektor[1])
            z.append(position_vektor[2])
            #energi=start_energi - energideponering

        elif energi<=0:
            break
            
        else:
            break
            # print('Partikel utanför sfär!')
    energideponering=start_energi-energi
    #plt.show()

    #
    return energideponering, x, y, z ,energideponering_list # , trajectory

#@jit(nopython=True)
def run_MC_alpha(iterationer, position_start_alpha, radie, max_antal_steg):
    
    energideponering_summa = 0
    x_lista=[]
    y_lista=[]
    z_lista=[]
    dos_lista=[]
    utanför = 0
    

    if position_start_alpha == position_start_alpha_skal:
        start_energi=energi_start(At211_energi,At211_sannolikhet)
        for i in range(iterationer):
            theta, phi = riktning_alpha()

            if not pi / 2 < phi < 3 * pi / 2:
                # print('Utanför')
                utanför += 1
                energideponering = 0
            else:
                start_position = position_start_alpha(radie, phi, theta)
                _,Steglängd = Stopping_power_och_steglängd(start_energi) #steglängd_alpha(start_position, df_stopping_power)
                steglängd=Steglängd*10**(-2) #Från cm till m
                energideponering, x,y,z, dos = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                       max_antal_steg)
                x_lista+=x
                y_lista+=y
                z_lista+=z
                dos_lista+=dos
    else:
        for i in range(iterationer):
            start_energi=energi_start(At211_energi,At211_sannolikhet)
            theta, phi = riktning_alpha()
            start_position = position_start_alpha(radie, phi, theta)
            _,Steglängd = Stopping_power_och_steglängd(start_energi) #steglängd_alpha(start_position, df_stopping_power)
            steglängd=Steglängd*10**(-2) #Från cm till m
            energideponering, x,y,z, dos  = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                   max_antal_steg)
            x_lista+=x
            y_lista+=y
            z_lista+=z
            dos_lista+=dos
        
        energideponering_summa += energideponering
    #Plotta en figur som visar energideponeringen i hela tumören
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x_lista,y_lista,z_lista,c=dos_lista,cmap='plasma') 

    print('Dos', dos_lista)
    print('antal utanför: ', utanför)
    print('total energideponering: ', energideponering_summa*10**6 )
    print(f'\nEnergideponering per partikel: {energideponering_summa*10**6 / (iterationer)} eV / partikel')#omvandlar från MeV till eV

    #Visa figur
    #plt.colorbar()
    plt.show()
    return energideponering_summa



if __name__ == "__main__":
    iterationer = 10 ** 2
    dummy_iterationer = 10**1
    max_antal_steg = 10**3
    radie_sfär = 300 * 10 ** (-6)
    

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, position_start_alpha_skal, radie_sfär,
                                        max_antal_steg)
    

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_alpha(iterationer, position_start_alpha_skal, radie_sfär, max_antal_steg)

    end_time(start)

    radie_sfär = 1 * 10 ** (-3)
    iterationer = 10 ** 2
    dummy_iterationer = 10**1
    max_antal_steg = 10**3
    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_alpha(dummy_iterationer, position_start_alpha_innanför, radie_sfär,
                                                 max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_alpha(iterationer, position_start_alpha_innanför, radie_sfär,
                                                 max_antal_steg)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nRESULTAT\n----------------------------------------------------------------------\n')

    print(f'\nSkal: Energideponering per partikel: {energideponering_tot_skal / iterationer:.2f} eV / partikel')
    print(f'Innanför: Energideponering per partikel: {energideponering_tot_innanför / iterationer:.2f} eV / partikel')


#Får samma siffror 0 MeV från skalet och 5.860718166904508 MeV i tumören
#Stämmer inte , får samma oavsett antal iterationer