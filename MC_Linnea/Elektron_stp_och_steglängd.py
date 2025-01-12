from imports import *
def Stopping_power_och_steglängd_elektron(Elektron_energi,rho_medium,stopping_power_data): #i MeV för elektron energin
    #Öppnar text filen och läser in varje koloumn
    """
    Data_energi=np.loadtxt('MC_Linnea/Elekt_stp_range_data')[:,0]
    Data_stoppingpower=np.loadtxt('MC_Linnea/Elekt_stp_range_data')[:,1]
    Data_range=(np.loadtxt('MC_Linnea/Elekt_stp_range_data')[:,2])
    """

    #Läser in värderna från data:

    # Omvandlar från MeV till eV
    energi_MeV_list = stopping_power_data[:, 0]
    energi_list = np.array([(lambda x: x * 10**6)(x) for x in energi_MeV_list])

    #från MeV cm^2/g till eV m^2/kg
    stopping_power_list = stopping_power_data[:, 1]
    STP_list =np.array([(lambda x: x * 10**5)(x) for x in stopping_power_list])

    #från g/cm^2 till kg/m^2
    CSDA_g_per_cm2_list=stopping_power_data[:, 2]
    CSDA_list=np.array([(lambda x: x * 10**1)(x) for x in CSDA_g_per_cm2_list])

    #Tar index för närmaste energi på alfapartikeln
    diff = np.abs( energi_list- Elektron_energi)
    closest_indices = np.argsort(diff)[:2]
    
    #Närmsta värdena:
    energi_close = energi_list[closest_indices]
    Stopping_power_close= stopping_power_list[closest_indices]
    CSDA_close=CSDA_list[closest_indices]

    
    #Linjär interpolera och få fram stoppingpower och slumpmässig steglängd
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        STP= Stopping_power_close[0]*rho_medium
        Steglängd=CSDA_close[0]/rho_medium
        Tau=Steglängd*np.random.random()
    else:
        STP = (Stopping_power_close[0] + (Elektron_energi - energi_close[0]) * (
                    Stopping_power_close[1] - Stopping_power_close[0]) / (
                                    energi_close[1] - energi_close[0]))*rho_medium
        
        Steglängd=(CSDA_close[0] + (Elektron_energi - energi_close[0]) * (
                    CSDA_close[1] - CSDA_close[0]) / (
                                    energi_close[1] - energi_close[0]))/rho_medium
        Tau=Steglängd*np.random.random()
    return STP, Steglängd, Tau #stp i eV/m och steglängd,Tau i m

