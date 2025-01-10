from imports import *
#Riktning är uniform för alfa, och kommer gå i samm riktning utan att ändra på vinkeln efter varje kollision
#from Martin.upg1_sampla_riktning_och_steg_start import riktning_uniform # Går inte att få upp VARFÖÖÖÖR???!!!! AHHHHHHHHHHHHHH
#from Martin.upg1_sampla_energi_start import energi_start

#Energi=energi_start(At211_energi,At211_sannolikhet) #Samplar på samma sätt alfa energin som fotonerna, är i MeV
from Elektron_stp_och_steglängd import Stopping_power_och_steglängd_elektron
from Elektron_energi import Elektron_startenergi
from Elektron_polarvinkel import Elektron_riktning


#Visa funktioner går att sätta jit men inte alla
@jit(nopython=True)
def riktning_elektron():
    theta = np.arccos(-1 + 2 * np.random.rand())
    phi = 2 * pi * np.random.rand()
    return theta, phi

@jit(nopython=True)
def förflyttning(position_vektor, steg_vektor):
    position_vektor += steg_vektor
    return position_vektor


#@jit(nopython=True)
def energiförlust_elektron(energi, steg):
    # Implementera stopping power
    # print('WIP')
    STP,_,_=Stopping_power_och_steglängd_elektron(energi)
    energiförlust=STP*steg #i MeV
    #energiförlust = energi * 0.1
    energi -= energiförlust

    if energi <= 0:
        energi = 0

    return energi


@jit(nopython=True)
def position_start_elektron_innanför(radie_sfär, phi, theta):
    r = radie_sfär * np.random.rand()

    x = r * np.sin(theta) * np.cos(phi)
    # y = r * np.sin(theta) * np.sin(phi)
    # z = r * np.cos(theta)
    y = 0
    z = 0

    position_vektor = np.array([x, y, z])
    return position_vektor


@jit(nopython=True)
def position_start_elektron_skal(radie_sfär, phi, theta):
    radie_elektron=1*10**(-21) #test för tillfället, sök upp elektronens radie 
    r = radie_sfär - 0.5 * radie_elektron  # För att inte endast theta = pi ska ge utslag

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

    # Under energideponering efter varje steg och ändrar riktningen
    while True:
        energideponering=0
        #Ändrar på riktningen efter varje steg
        riktning = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.cos(phi), np.cos(theta)])
        riktning /= np.linalg.norm(riktning)

        #Stegstorleken när den ändrar riktning
        _,Steg,Tau=Stopping_power_och_steglängd_elektron(energi-energideponering)
        steg_storlek = (Steg-Tau)*10**(-2) #omvandlar cm till m

        steg_vektor = riktning * steg_storlek

        position_vektor += steg_vektor
        energi = energiförlust_elektron(energi, steg_storlek)
        
        
        if np.dot(position_vektor, position_vektor) <= radie:
            innanför = True
            # trajectory.append(tuple(position_vektor))
            #print(f'Energideponering i position ', position_vektor)
            energideponering += start_energi - energi
            theta=Elektron_riktning(start_energi-energideponering)
            phi=np.random.random()*2*pi

        elif energi-energideponering<=0:
            #Förlorat all sin energi
            break
            
        else:
            break
            # print('Partikel utanför sfär!')

    

    return energideponering  # , trajectory

#@jit(nopython=True)
def run_MC_elektron(iterationer, position_start_eletron, radie, max_antal_steg):
    
    energideponering_summa = 0
    utanför = 0
    start_energi=Elektron_startenergi()

    if position_start_eletron == position_start_elektron_skal:

        for i in range(iterationer):
            theta, phi = riktning_elektron()

            if not pi / 2 < phi < 3 * pi / 2:
                # print('Utanför')
                utanför += 1
                energideponering = 0
            else:
                start_position = position_start_eletron(radie, phi, theta)
                _,_,Steglängd = Stopping_power_och_steglängd_elektron(start_energi) #steglängd_alpha(start_position, df_stopping_power)
                steglängd=Steglängd*10**(-2)
                energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                       max_antal_steg)

    else:
        for i in range(iterationer):
            theta, phi = riktning_elektron()
            start_position = position_start_eletron(radie, phi, theta)
            _,_,steglängd = Stopping_power_och_steglängd_elektron(start_energi) #steglängd_alpha(start_position, df_stopping_power)
            energideponering = laddad_partikel_väg(start_energi, start_position, phi, theta, steglängd, radie,
                                                   max_antal_steg)

        energideponering_summa += energideponering

    print('antal utanför: ', utanför)
    print('total energideponering: ', energideponering_summa)
    print(f'\nEnergideponering per partikel: {energideponering_summa / (iterationer*10**6):.2f} eV / partikel')
    return energideponering_summa



if __name__ == "__main__":
    iterationer = 10 ** 2
    dummy_iterationer = 10**2
    max_antal_steg = 10**3
    radie_sfär = 300 * 10 ** (-6)
    

    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_elektron(dummy_iterationer, position_start_elektron_skal, radie_sfär,
                                        max_antal_steg)
    

    start = time.time()

    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_skal = run_MC_elektron(iterationer, position_start_elektron_skal, radie_sfär, max_antal_steg)

    end_time(start)

    radie_sfär = 1 * 10 ** (-3)
    iterationer = 10 ** 2
    dummy_iterationer = 10**1
    max_antal_steg = 10**3
    print(
        '\n----------------------------------------------------------------------\nDUMMY\n----------------------------------------------------------------------\n')

    _ = run_MC_elektron(dummy_iterationer, position_start_elektron_innanför, radie_sfär,
                                                 max_antal_steg)

    start = time.time()
    print(
        '\n----------------------------------------------------------------------\nRIKTIG\n----------------------------------------------------------------------\n')
    energideponering_tot_innanför = run_MC_elektron(iterationer, position_start_elektron_innanför, radie_sfär,
                                                 max_antal_steg)

    end_time(start)

    print(
        '\n----------------------------------------------------------------------\nRESULTAT\n----------------------------------------------------------------------\n')

    print(f'\nSkal: Energideponering per partikel: {energideponering_tot_skal / iterationer:.2f} eV / partikel')
    print(f'Innanför: Energideponering per partikel: {energideponering_tot_innanför / iterationer:.2f} eV / partikel')


#Får samma siffror 0 MeV från skalet och 5.860718166904508 MeV i tumören
#Stämmer inte , får samma oavsett antal iterationer