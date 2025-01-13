from imports import *
from upg2_energi_efter_förlust_elektron import energi_efter_energiförlust
from upg2_steg_transformation import ny_steg_transformera_koordinatsystem_3d
from upg2_Elektron_stp_och_steglängd import Stopping_power_och_steglängd_elektron
from upg2_Elektron_polarvinkel import Elektron_theata_ny

def laddad_partikel_väg_elektron(energi_start, position_start, phi, theta, steglängd, radie_sfär, rho_medium,
                        stopping_power_data,
                        max_antal_steg):
    """
    Funktion som följer alfapartikeln allteftersom den växelverkar i ett medium.
    :param radie_sfär: Radien av sfären för fördelningen.
    :param max_antal_steg: Maximalt antal steg som steglängden ska delas upp i.
    :return: Energideponeringen innanför sfären.
    """

    position_vektor = position_start
    energi = energi_start
    steg_tagna = 0
    x,y,z,dos=[],[],[],[]                                
    stopping_power_data = np.loadtxt('MC_Linnea/Elekt_stp_range_data')
    scatter_power_data=np.loadtxt('MC_Linnea/Scatterpower_vatten_data')
    # Under tiden som partikeln fortfarande inte tagit hela sitt steg.
    while steg_tagna < max_antal_steg and energi > 0:
        if np.dot(position_vektor, position_vektor) > radie_sfär**2:
            break
        riktning = np.array(
        [np.sin(theta) * np.cos(phi)
            , np.sin(theta) * np.cos(phi)
            , np.cos(theta)])

        riktning /= np.linalg.norm(riktning)
        position_vektor += riktning * steglängd

        #Nya värden på vinklarna
        phi_ny=np.random.random()*2*pi
        theta_ny=Elektron_theata_ny(energi,scatter_power_data,rho_vatten)

        #stegstorlek totalt blir steglängd+cos(theta)*(s-steglängd) där s är CSDA
        _,s,_=Stopping_power_och_steglängd_elektron(energi,rho_medium,stopping_power_data)

        #Ändrar på positionsvektor efter att transformations matrisen
        position_vektor+=ny_steg_transformera_koordinatsystem_3d(steglängd,theta,phi,s-steglängd,phi_ny,theta_ny)
        
        #Plottvärderna för att se dosfördelningen, men får bara ut startpositionen inte alla andra delsteg...
        dos.append(energi-energi_efter_energiförlust(energi, steglängd, rho_medium, stopping_power_data))
        x.append(position_vektor[0])
        y.append(position_vektor[1])
        z.append(position_vektor[2])

        #Förlorar energi och får en ny efter endast steglängden
        energi = energi_efter_energiförlust(energi, steglängd, rho_medium, stopping_power_data)

        #Ändrar på vinklarna och värdet på steglängden innan nästa vinkeländring
        phi=phi_ny
        theta=theta_ny
        steglängd=s-steglängd

        #Tar ett steg
        steg_tagna += 1
        
        # if np.dot(position_vektor, position_vektor) <= radie_sfär:
        #print(f'Energideponering i position ', position_vektor)
        # else:
        #     break
        print

    energideponering = energi_start - energi
    #print('Doslista', np.sum(dos))
    #print(f'energideponering: {energideponering} eV')
    return energideponering , x,y,z, dos
