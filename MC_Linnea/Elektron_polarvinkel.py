from imports import *


def Elektron_theata_ny(Elektron_energi): #Energi i elektron i tabell är i eV
    R=np.random.random() #Slumpmässig tal mellan 0-1

    #Anta att tumören är i en vävnad alltså se på vatten för tvärsnitten

    #Läser in excel filen och tar ut koloumn med energi och scatter power
    file_scatter=pd.read_excel(r'MC_Linnea/Monte Carlo, test.xlsx', sheet_name='Scattering power') #Behövde skriva av en egen Excel fil från ICRU 35
    Energi_data=file_scatter['E/MeV'].to_list()
    Mass_Scatterpower_data=file_scatter['Scatter power Water '].to_list()  #T/rho för vatten
    
    #Hittar närmast energi som liknar elektron energin och ta fram scatter power:
    rho=0.998 #densitet på vatten i g/cm^3

    #Tar index för närmaste energi på elektronen
    energi_list = np.array(Energi_data) #Gör om listorna till arrays
    Mass_Scatterpower_list=np.array(Mass_Scatterpower_data) #Gör om listorna till arrays
    diff = np.abs(energi_list - Elektron_energi)
    closest_indices = np.argsort(diff)[:2]


    #Få scattering power nära energierna och energivärderna närmast

    T_close=Mass_Scatterpower_list[closest_indices]*rho
    energi_close = energi_list[closest_indices]



    #Linjär interpolera och få fram theata_s
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        theata_s= T_close[0]

    else:
        theata_s = T_close[0] + (Elektron_energi - energi_close[0]) * (
                    T_close[1] - T_close[0]) / (
                                    energi_close[1] - energi_close[0])
 

    return np.sqrt(-theata_s*np.log(1-R))


    
