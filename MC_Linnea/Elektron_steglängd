from imports import *
#from Elektron_polarvinkel import Elektron_riktning #Ta reda på varför det inte fungerar, är i samma mapp 

def Steglängd_elektron(Elektron_energi): #Energi i MeV    
    #Läser in excelfil med data 
    Data=pd.read_excel(r'given_data/Tvärsnittstabell_Elektroner.xlsx')

    #Tar ut varje koloumn med data
    Energi_data=Data['Energy (eV)'].to_list() #i eV gör om till MeV eller tvärtom
    Range_data=Data['Range (g/cm^2)'].to_list()

    #Omvandlar till enhet i MeV
    Energi_data=[Energi_data[i]/10**6 for i in range(len(Energi_data))] 

    #Gör till listor
    Energi_list=np.array(Energi_data)
    Range_list=np.array(Range_data)

    #Närmaste index till elektronens energi
    diff = np.abs(Energi_list - Elektron_energi)
    closest_indices = np.argsort(diff)[:2]

    #Närmaste elektronens energi och motsvarande range:
    Energi_close=Energi_list[closest_indices]
    Range_close=Range_list[closest_indices]

    #Linjär interpolera och få fram theata_s
    if Energi_close[1] - Energi_close[0] < 10 ** (-15):
        Steglängd= Range_close[0]
        Tau=np.random.random()*Steglängd
        

    else:
        Steglängd = Range_close[0] + (Elektron_energi - Energi_close[0]) * (
                    Range_close[1] - Range_close[0]) / (
                                    Energi_close[1] - Energi_close[0])
        Tau=np.random.random()*Steglängd
 
    return Steglängd-Tau
    

