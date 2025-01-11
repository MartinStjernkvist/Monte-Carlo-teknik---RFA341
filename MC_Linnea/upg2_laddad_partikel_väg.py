from imports import *
from upg2_energi_efter_förlust import energi_efter_energiförlust


def laddad_partikel_väg(energi_start, position_start, phi, theta, steglängd, radie_sfär, rho_medium,
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

    steg_storlek = steglängd / max_antal_steg

    riktning = np.array(
        [np.sin(theta) * np.cos(phi)
            , np.sin(theta) * np.cos(phi)
            , np.cos(theta)])

    riktning /= np.linalg.norm(riktning)
    steg_vektor = riktning * steg_storlek

    steg_tagna = 0
    x,y,z,dos=[],[],[],[]
    # Under tiden som partikeln fortfarande inte tagit hela sitt steg.
    while steg_tagna < max_antal_steg and energi > 0:
        #Plottvärderna för att se dosfördelningen
        x.append(position_vektor[0])
        y.append(position_vektor[1])
        z.append(position_vektor[2])
        dos.append(energi-energi_efter_energiförlust(energi, steg_storlek, rho_medium, stopping_power_data))
        
        steg_tagna += 1
        position_vektor += steg_vektor

        energi = energi_efter_energiförlust(energi, steg_storlek, rho_medium, stopping_power_data)

        
        # if np.dot(position_vektor, position_vektor) <= radie_sfär:
        #     print(f'Energideponering i position ', position_vektor)
        # else:
        #     break
        if np.dot(position_vektor, position_vektor) > radie_sfär:
            break
    energideponering = energi_start - energi
    #print('Doslista', np.sum(dos))
    print(f'energideponering: {energideponering} eV')
    return energideponering#, x,y,z, dos
