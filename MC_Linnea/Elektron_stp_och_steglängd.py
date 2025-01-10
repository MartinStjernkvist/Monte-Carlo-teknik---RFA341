from imports import *
def Stopping_power_och_steglängd_elektron(Elektron_energi): #i MeV för elektron energin
    #Öppnar text filen och läser in varje koloumn
    Data_energi=np.loadtxt('MC_Linnea/Elekt_stp_range_data')[:,0]
    Data_stoppingpower=np.loadtxt('MC_Linnea/Elekt_stp_range_data')[:,1]
    Data_range=(np.loadtxt('MC_Linnea/Elekt_stp_range_data')[:,2])

    #Tar index för närmaste energi på alfapartikeln
    diff = np.abs(Data_energi - Elektron_energi)
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
        Tau=Steglängd*np.random.random()
    else:
        STP = (Stopping_power_close[0] + (Elektron_energi - energi_close[0]) * (
                    Stopping_power_close[1] - Stopping_power_close[0]) / (
                                    energi_close[1] - energi_close[0]))*rho
        
        Steglängd=(Range_close[0] + (Elektron_energi - energi_close[0]) * (
                    Range_close[1] - Range_close[0]) / (
                                    energi_close[1] - energi_close[0]))/rho
        Tau=Steglängd*np.random.random()
    return STP, Steglängd, Tau #stp i MeV/cm och steglängd i cm

