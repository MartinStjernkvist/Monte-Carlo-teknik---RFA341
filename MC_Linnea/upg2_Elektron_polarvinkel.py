from imports import *


def Elektron_theata_ny(Elektron_energi, Scatterpower_data, rho_medium): #Energi i elektron i tabell är i eV
    R=np.random.random() #Slumpmässig tal mellan 0-1

    #Anta att tumören är i en vävnad alltså se på vatten för tvärsnitten

    """
    #Läser in excel filen och tar ut koloumn med energi och scatter power
    file_scatter=pd.read_excel(r'MC_Linnea/Monte Carlo, test.xlsx', sheet_name='Scattering power') #Behövde skriva av en egen Excel fil från ICRU 35
    Energi_data=file_scatter['E/MeV'].to_list()
    Mass_Scatterpower_data=file_scatter['Scatter power Water '].to_list()  #T/rho för vatten
    """
  
    # Omvandlar från MeV till eV och till en np.array
    energi_MeV_list = Scatterpower_data[:,0]
    energi_list = np.array([(lambda x: x * 10**6)(x) for x in energi_MeV_list])

    Mass_Scatterpower_data=Scatterpower_data[:,1]
    Mass_Scatterpower_list=np.array(Mass_Scatterpower_data) #Gör om listorna till arrays

    #Hittar närmast energi som liknar elektron energin och ta fram scatter power:

    #Tar index för närmaste energi på elektronen
    diff = np.abs(energi_list - Elektron_energi)
    closest_indices = np.argsort(diff)[:2]

    #Få scattering power nära energierna och energivärderna närmast
    T_close=Mass_Scatterpower_list[closest_indices]*rho_medium*10**(-3) #kg/m^3 till g/cm^3
    energi_close = energi_list[closest_indices]

    #Linjär interpolera och få fram theata_s
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        theata_s= T_close[0]

    else:
        theata_s = T_close[0] + (Elektron_energi - energi_close[0]) * (
                    T_close[1] - T_close[0]) / (
                                    energi_close[1] - energi_close[0])
 

    return np.sqrt(-theata_s*np.log(1-R))


    
