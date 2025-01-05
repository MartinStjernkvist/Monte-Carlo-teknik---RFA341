from imports import *



class attenueringsdata:

    def __init__(self, voxelvärde, energi, df_attenueringsdata, df_anatomidefinitioner):

        # Läs in data från argumenten till klassen. Excelfilerna bearbetas med paketet pandas.
        self.df_attenueringsdata = df_attenueringsdata
        self.anatomidefinitioner = df_anatomidefinitioner

        self.voxelvärde = voxelvärde
        self.energi = energi
        self.energi_list = self.df_attenueringsdata['E'].to_list()

    def voxelvärde_till_material(self):
        """
        Funktion som omvandlar voxelvärdet i fantommatrisen till ett material.
        Voxelvärdet motsvarar viss sorts vävnad, enligt filen Anatomidefinitioner.xlsx.
        Därför behövs en koppling mellan vävnaderna i anatomidefinitioner och materialen i filen Attenueringsdata.xlsx.
        :return: 1) materialnamnet som motsvarar ett visst voxelvärde och 2) en namnlista med material som ingår i fantomen.
        """
        # print(vävnad)

        # luft = 0
        # bihåla? antag vatten
        # matstrupe? antag vatten
        # svalg? antag muskler
        # kortikalt ben? antag röd benmärg
        # benmärg? antingen röd eller yellow, antag röd
        # prostata? antar bladder
        # skull? antag benmärg, antag röd
        # bronker? antag vatten

        # Skapa listor för materialen, där ett material motsvarar voxelvärden (vävnad enligt Anatomidefinitioner).
        water = [3, 8, 12, 37]
        muscle = [6, 13]
        lung = [11]
        dry_spine = [27, 28]
        dry_rib = [25]
        # adipose = [] # Hittade ingenting som passade in här
        blood = [2]
        heart = [1]
        kidney = [17, 18, 19, 20, 21, 22, 23]
        liver = [9]
        lymph = [36, 41]
        pancreas = [16]
        intestine = [14, 15, 32, 33]
        # skull = [] # Hittade ingenting som passade in här
        cartilage = [34]
        brain = [7]
        spleen = [24]
        air = [0, 35, 38]
        breast_mammary = [5]
        skin = [4]
        eye_lens = [42, 43]
        # ovary = [] # Hittade ingenting som passade in här
        red_marrow = [26, 29]
        yellow_marrow = [44]
        # testis = [] # Hittade ingenting som passade in här
        thyroid = [39, 40]
        bladder = [10, 30, 31]

        # Sätt ihop individuella materiallistor till en stor lista.
        big_list = [water, muscle, lung, dry_spine, dry_rib, blood, heart, kidney, liver, lymph, pancreas, intestine,
                    cartilage, brain, spleen, air, breast_mammary, skin, eye_lens, red_marrow, yellow_marrow, thyroid,
                    bladder]

        # Lista med namn, för att kunna indexera i attenueringsdatan.
        big_list_names = ['water', 'muscle', 'lung', 'dry_spine', 'dry_rib', 'blood', 'heart', 'kidney', 'liver',
                          'lymph', 'pancreas', 'intestine',
                          'cartilage', 'brain', 'spleen', 'air', 'breast_mammary', 'skin', 'eye_lens', 'red_marrow',
                          'yellow_marrow', 'thyroid',
                          'bladder']

        # Filtrerar en materiallista i taget i den större, sammansatta listan.
        # Om voxelvärdet matchar ett av värdena i en materiallista erhålls materialnamnet.
        for i in range(len(big_list)):
            sublist = big_list[i]
            if any(self.voxelvärde == värde for värde in sublist):
                material = big_list_names[i]
                break

        return material, big_list_names

    def mu(self):
        """
        Funktion som tar fram attenueringskoefficienten utifrån energin och voxelvärdet i en position.
        :return: attenueringskoefficient
        """
        # Energilistan (första kolumnen i Attenueringsdata.xlsx) läses in.
        # De två närmaste indexen, för energierna i listan som är närmast den ingående energin för funktionen tas fram.
        energi_list = np.array(self.energi_list)
        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]

        # Kalla på ovanstående funktion, för att omvandla voxelvärde till material.
        material, _ = self.voxelvärde_till_material()
        mu_list = np.array(self.df_attenueringsdata[material].to_list())

        # Närmaste energierna och attenueringskoefficienterna.
        energi_close = energi_list[closest_indices]
        mu_close = mu_list[closest_indices]

        # Om närliggande attenueringskoefficienter är lika stora -> attenueringskoefficient = första värdet.
        if mu_close[1] - mu_close[0] < 10 ** (-15):
            mu_target = mu_close[0]

        # Om närliggande attenueringskoefficienter inte är lika stora -> linjärinterpolera fram attenueringskoefficienten.
        else:
            mu_target = mu_close[0] + (self.energi - energi_close[0]) * (mu_close[1] - mu_close[0]) / (
                    energi_close[1] - energi_close[0])

        # print(f'mu target {mu_target}')

        return mu_target


    def mu_max(self):
        """
        Funktion som tar fram den maximala attenueringskoefficienten i fantomen för en viss energi.
        :return: maximal attenueringskoefficient
        """
        # Energilistan (första kolumnen i Attenueringsdata.xlsx) läses in.
        # De två närmaste indexen, för energierna i listan som är närmast den ingående energin för funktionen tas fram.
        energi_list = np.array(self.energi_list)
        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]

        # Kalla på ovanstående funktion, för att erhålla namnlistan med material som ingår i fantomen.
        _, big_list_names = self.voxelvärde_till_material()

        # Skapa en tom vektor, som sedan ska fyllas.
        mu_array = np.zeros(len(big_list_names), dtype=object)

        # För varje material i namnlistan tas värden på attenueringskoefficienten för samma "närmaste index" som för energin.
        # En vektor med listor skapas.
        for i in range(len(big_list_names)):
            mu_array[i] = np.array(self.df_attenueringsdata[big_list_names[i]][closest_indices].to_list())

        # Sortera bland vektorn med listorna, för att ta fram materialet (indexet) som har störst attenueringskoefficient.
        mu_max_close_values = [np.max(arr) for arr in mu_array]
        mu_max_close = max(mu_max_close_values)
        mu_max_close_index = mu_max_close_values.index(mu_max_close)

        # Ta ut listan med de två attenueringskoefficienterna, för materialet med högst attenueringskoefficient.
        mu_max_close = mu_array[mu_max_close_index]
        # Även energierna som ligger närmast energin som matas in i funktionen.
        energi_close = energi_list[closest_indices]

        # Om närliggande attenueringskoefficienter är lika stora -> attenueringskoefficient = första värdet.
        if mu_max_close[1] - mu_max_close[0] < 10 ** (-15):
            mu_max = mu_max_close[0]

        # Om närliggande attenueringskoefficienter inte är lika stora -> linjärinterpolera fram attenueringskoefficienten.
        else:
            mu_max = mu_max_close[0] + (self.energi - energi_close[0]) * (mu_max_close[1] - mu_max_close[0]) / (
                    energi_close[1] - energi_close[0])

        return mu_max



if __name__ == "__main__":
    start = time.time()

    df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)
    # print(df_attenueringsdata.columns)
    # print(df_attenueringsdata['muscle'])
    #
    # energi_list = df_attenueringsdata['E'].to_list()
    # print(energi_list)

    df_anatomidefinitioner = pd.read_excel(anatomidefinitioner_file, index_col=None)
    # print(df_anatomidefinitioner.columns)
    #
    # print(df_anatomidefinitioner)
    # print(df_anatomidefinitioner['Unnamed: 1'])

    voxelvärde = 1
    energi = 10000

    instans = attenueringsdata(voxelvärde, energi, df_attenueringsdata, df_anatomidefinitioner)
    mu = instans.mu()
    print(mu)

    end_time(start)
