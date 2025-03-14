from imports import *
def Stopping_power_och_steglängd(Alfa_energi): #i MeV för alfa energin
    #Öppnar text filen och läser in varje koloumn
    Data_energi=np.loadtxt('MC_Linnea/Stoppingpower_data_alfa')[:,0]
    Data_stoppingpower=np.loadtxt('MC_Linnea/Stoppingpower_data_alfa')[:,1]
    Data_range=(np.loadtxt('MC_Linnea/Stoppingpower_data_alfa')[:,2])

    #Tar index för närmaste energi på alfapartikeln
    diff = np.abs(Data_energi - Alfa_energi)
    closest_indices = np.argsort(diff)[:2]
    
    #Närmsta värdena:
    energi_close = Data_energi[closest_indices]
    Stopping_power_close=Data_stoppingpower[closest_indices]
    Range_close=Data_range[closest_indices]

    rho= 0.998 #vattens densitet i g/cm^3
    #Linjär interpolera och få fram stoppingpower och slumpmässig steglängd
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        STP= Stopping_power_close[0]*rho
        Steglängd=Range_close[0]/rho

    else:
        STP = (Stopping_power_close[0] + (Alfa_energi - energi_close[0]) * (
                    Stopping_power_close[1] - Stopping_power_close[0]) / (
                                    energi_close[1] - energi_close[0]))*rho
        
        Steglängd=(Range_close[0] + (Alfa_energi - energi_close[0]) * (
                    Range_close[1] - Range_close[0]) / (
                                    energi_close[1] - energi_close[0]))/rho
        
    return STP*10**8, Steglängd*10**(-2) #stp i eV/m och steglängd i m

