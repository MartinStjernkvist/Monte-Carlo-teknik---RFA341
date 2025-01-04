from imports import *


class växelverkan:

    def __init__(self, energi, df_tvärsnitt):
        self.df = df_tvärsnitt

        self.energi = energi
        self.energi_list = self.df['Energy (eV)'].to_list()
        self.foto_list = self.df['Photoelectric  (cm^2)'].to_list()
        self.compton_list = self.df['Compton (cm^2)'].to_list()
        self.rayleigh_list = self.df['Rayleigh (cm^2)'].to_list()

    def find_foto_tvärsnitt(self):
        energi_list = np.array(self.energi_list)
        foto_list = np.array(self.foto_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        foto_close = foto_list[closest_indices]
        energi_close = energi_list[closest_indices]

        if energi_close[1] - energi_close[0] < 10 ** (-15):
            foto_target = foto_close[0]

        else:
            # linjär interpolering funktion
            foto_target = foto_close[0] + (self.energi - energi_close[0]) * (foto_close[1] - foto_close[0]) / (
                    energi_close[1] - energi_close[0])
        return foto_target

    def find_compton_tvärsnitt(self):
        energi_list = np.array(self.energi_list)
        compton_list = np.array(self.compton_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        compton_close = compton_list[closest_indices]
        energi_close = energi_list[closest_indices]

        if energi_close[1] - energi_close[0] < 10 ** (-15):
            compton_target = compton_close[0]

        else:
            compton_target = compton_close[0] + (self.energi - energi_close[0]) * (compton_close[1] - compton_close[0]) / (
                    energi_close[1] - energi_close[0])

        return compton_target


    def find_rayleigh_tvärsnitt(self):
        energi_list = np.array(self.energi_list)
        rayleigh_list = np.array(self.rayleigh_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        rayleigh_close = rayleigh_list[closest_indices]
        energi_close = energi_list[closest_indices]

        if energi_close[1] - energi_close[0] < 10**(-15):
            rayleigh_target = rayleigh_close[0]

        else:
            rayleigh_target = rayleigh_close[0] + (self.energi - energi_close[0]) * (rayleigh_close[1] - rayleigh_close[0]) / (
                    energi_close[1] - energi_close[0])
        return rayleigh_target

    # @timer
    def bestäm_växelverkan(self):
        foto_target = self.find_foto_tvärsnitt()
        compton_target = self.find_compton_tvärsnitt()
        rayleigh_target = self.find_rayleigh_tvärsnitt()

        tvärsnitt_lista = [foto_target, compton_target, rayleigh_target]  # cm^2

        tvärsnitt_lista_norm = np.cumsum(tvärsnitt_lista) / np.sum(tvärsnitt_lista)

        slump_tal = np.random.rand()
        if slump_tal <= tvärsnitt_lista_norm[0]:
            text = 'foto'
        elif slump_tal <= tvärsnitt_lista_norm[1]:
            text = 'compton'
        else:
            text = 'rayleigh'
        return text



class växelverkan_slimmad:

    def __init__(self, energi, df_tvärsnitt):
        self.df = df_tvärsnitt

        self.energi = energi
        self.energi_list = self.df['Energy (eV)'].to_list()
        self.foto_list = self.df['Photoelectric  (cm^2)'].to_list()
        self.compton_list = self.df['Compton (cm^2)'].to_list()
        self.rayleigh_list = self.df['Rayleigh (cm^2)'].to_list()

    def find_foto_tvärsnitt_slimmad(self):
        energi_list = np.array(self.energi_list)
        foto_list = np.array(self.foto_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        foto_close = foto_list[closest_indices]

        foto_target = foto_close[0]

        return foto_target

    def find_compton_tvärsnitt_slimmad(self):
        energi_list = np.array(self.energi_list)
        compton_list = np.array(self.compton_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        compton_close = compton_list[closest_indices]

        compton_target = compton_close[0]

        return compton_target

    def find_rayleigh_tvärsnitt_slimmad(self):
        energi_list = np.array(self.energi_list)
        rayleigh_list = np.array(self.rayleigh_list)

        diff = np.abs(energi_list - self.energi)
        closest_indices = np.argsort(diff)[:2]
        rayleigh_close = rayleigh_list[closest_indices]

        rayleigh_target = rayleigh_close[0]

        return rayleigh_target


    def bestäm_växelverkan_slimmad(self):
        foto_target = self.find_foto_tvärsnitt_slimmad()
        compton_target = self.find_compton_tvärsnitt_slimmad()
        rayleigh_target = self.find_rayleigh_tvärsnitt_slimmad()

        tvärsnitt_lista = [foto_target, compton_target, rayleigh_target]  # cm^2

        tvärsnitt_lista_norm = np.cumsum(tvärsnitt_lista) / np.sum(tvärsnitt_lista)

        slump_tal = np.random.rand()
        if slump_tal <= tvärsnitt_lista_norm[0]:
            text = 'foto'
        elif slump_tal <= tvärsnitt_lista_norm[1]:
            text = 'compton'
        else:
            text = 'rayleigh'
        return text

#   ----------------------------------------------------------------------
#   FOTO OCH COMPTON
#   ----------------------------------------------------------------------


if __name__ == "__main__":
    start = time.time()

    #   ----------------------------------------------------------------------
    #   INPUT DATA
    #   ----------------------------------------------------------------------
    df_tvärsnitt = pd.read_excel(tvärsnitt_file, index_col=None)
    # print(df.columns)

    energi_list = df_tvärsnitt['Energy (eV)'].to_list()
    compton_list = df_tvärsnitt['Compton (cm^2)'].to_list()
    foto_list = df_tvärsnitt['Photoelectric  (cm^2)'].to_list()
    rayleigh_list = df_tvärsnitt['Rayleigh (cm^2)'].to_list()

    #   ----------------------------------------------------------------------
    #   INPUT ENERGI
    #   ----------------------------------------------------------------------
    energi = 10 ** 3

    #   ----------------------------------------------------------------------
    #   LÅT STÅ
    #   ----------------------------------------------------------------------

    instans = växelverkan(energi, df_tvärsnitt)

    foto_target = instans.find_foto_tvärsnitt()
    compton_target = instans.find_compton_tvärsnitt()
    rayleigh_target = instans.find_rayleigh_tvärsnitt()

    x_data = [energi_list, energi, energi_list, energi, energi_list, energi]
    y_data = [foto_list, foto_target, compton_list, compton_target, rayleigh_list, rayleigh_target]
    scatter = [2, 1, 2, 1, 2, 1]
    label_data = ['foto', 'foto', 'compton', 'compton', 'rayleigh', 'rayleigh']
    marker = ['o', 'X', 'o', 'X', 'o', 'X']
    color = ['blue', 'red', 'green', 'red', 'magenta', 'red']

    fig = plot_stuff(x_data, y_data, scatter, label_data,
                     marker, color, x_label='energi (eV)', y_label='tvärsnitt (cm^2)', title='foto och compton',
                     fig_size=(10, 10), symbol_size=100, font_size=30, alpha=1, line_width=2, x_lim=(10, 10**6),
                     y_lim=(10**(-30), 10**(-17)),
                     grid=True, x_scale='log', y_scale='log')

    fig.savefig('tvärsnitt', bbox_inches='tight')

    print(instans.bestäm_växelverkan())

    # iterationer = 1000
    # bingo = 0
    # for i in range(iterationer):
    #     instans = växelverkan(energi, tvärsnitt_file)
    #     if instans.bestäm_växelverkan() == 'foto':
    #         bingo += 1
    # print(bingo)

    end_time(start)

"""
# INSÅG INTE ATT VI HADE EN EXCELFIL MED TVÄRSNITT

class växelverkan:

    def sigma_compton(self, energi):
        # klein nishina, podgorsak
        epsilon = energi / E_e
        Z = 10
        if energi < 100_000:
            sigma_compton = Z * (8 / 3) * pi * r_e ** 2 * 1 / (1 + 2 * epsilon) ** 2 * (
                    1 + 2 * epsilon + (6 / 5) * epsilon ** 2 - (1 / 2) * epsilon ** 3 + (2 / 7) * epsilon ** 4 - (
                    6 / 35) * epsilon ** 5 + (8 / 105) * epsilon ** 6 + (4 / 105) * epsilon ** 7)
        else:
            sigma_compton = Z * 2 * pi * r_e ** 2 * ((1 + epsilon) / epsilon ** 2 * (
                    2 * (1 + epsilon) / (1 + 2 * epsilon) - np.log(1 + 2 * epsilon) / epsilon) + np.log(
                1 + 2 * epsilon) / (2 * epsilon) - (1 + 3 * epsilon) / (1 + 2 * epsilon) ** 2)

        return sigma_compton

    def sigma_foto(self, energi):
        # https://en.wikipedia.org/wiki/Gamma_ray_cross_section
        Z = 10  # räknar inte med Z-beroende, BYT TILL MEDELVÄRDE AV Z FÖR MÄNNISKOKROPP
        alpha = 1 / 137
        k = energi / E_e
        return 16 / 3 * np.sqrt(2) * pi * r_e ** 2 * alpha ** 4 * Z ** 5 / k ** 3.5
        # return 5 * 10**11 * Z**5 / energi**(3.5)
        # return 16 / 3 * np.sqrt(2) * pi * alpha ** 8 * a_0 ** 2 * (m_e * c ** 2) ** 3.5 * Z ** 5 / energi ** 3.5

    def bestäm_växelverkan(self, energi):
        # energierna på gamma från sönderfallen är för låga för att ge parbildning
        # -> antingen compton eller foto
        # möjligtvis även rayleigh
        sigma_foto = växelverkan.sigma_foto(self, energi)
        sigma_compton = växelverkan.sigma_compton(self, energi)

        tvärsnitt_lista = [sigma_foto, sigma_compton]  # keV
        # print(tvärsnitt_lista)

        tvärsnitt_lista_norm = np.zeros(len(tvärsnitt_lista))
        for i in range(len(tvärsnitt_lista)):
            if i == 0:
                tvärsnitt_lista_norm[i] = tvärsnitt_lista[i] / np.sum(tvärsnitt_lista)
            else:
                tvärsnitt_lista_norm[i] = tvärsnitt_lista[i] / np.sum(tvärsnitt_lista) + tvärsnitt_lista_norm[i - 1]
        # print(tvärsnitt_lista_norm)

        # OBS, kanske måste räkna med rayleigh spridning
        if np.random.rand() <= tvärsnitt_lista[0]:
            text = 'foto'
        else:
            text = 'compton'

        # print(text)
        return text

if __name__ == "__main__":
    start = 10_000
    stop = 500_000

    # x_data = np.linspace(start, stop)
    # y_data = växelverkan().sigma_foto(x_data)
    # scatter = 2
    # label_data = 'tvärsnitt foto'
    #
    # fig = plot_stuff(x_data, y_data, scatter, label_data,
    #                  marker='o', color='green', x_label='energi (eV)', y_label='sigma (barn)', title='tvärsnitt foto',
    #                  fig_size=(10, 10), symbol_size=50, font_size=30, alpha=1, line_width=5, x_lim=(0, 0), y_lim=(0, 0),
    #                  grid=True, x_scale='log', y_scale='log')
    #
    # fig.savefig('foto.png', bbox_inches='tight')


    x_data = [np.linspace(start, stop), np.linspace(start, stop)]
    y_data = [list(map(växelverkan().sigma_foto, x_data[0])), list(map(växelverkan().sigma_compton, x_data[1]))]
    scatter = [2, 2]
    label_data = ['foto', 'compton']
    color = ['blue', 'red']

    fig = plot_stuff(x_data, y_data, scatter, label_data,
                     marker='o', color=color, x_label='energi (eV)', y_label='sigma (barn)', title='tvärsnitt foto',
                     fig_size=(10, 10), symbol_size=50, font_size=30, alpha=1, line_width=5, x_lim=(start, stop),
                     y_lim=(10 ** (-1), 10 ** 5),
                     grid=True, x_scale='log', y_scale='log')

    fig.savefig('foto & compton.png', bbox_inches='tight')

    iterationer = 1000
    bingo = []
    for i in range(iterationer):
        instance = växelverkan().bestäm_växelverkan(100_000)
        if instance == 'foto':
            bingo.append(1)

    print(len(bingo))
"""
