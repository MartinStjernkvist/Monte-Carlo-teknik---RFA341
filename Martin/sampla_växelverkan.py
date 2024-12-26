from imports import *

tvärsnitt_file = '../given_data/Tvärsnittstabeller_Fotoner.xlsx'
df = pd.read_excel(tvärsnitt_file, index_col=None)
# print(df)
print(df.columns)
# print(df.at[5, 'Unnamed: 1'])
# print(df['Rayleigh (cm^2)'])

# first_column = df.iloc[:, 0]
# print(first_column)
energi_list = df['Energy (eV)'].to_list()
# print(energi_list)

compton_list = df['Compton (cm^2)'].to_list()
# print(compton_list)

foto_list = df['Photoelectric  (cm^2)'].to_list()
# print(foto_list)

initial_guess = [-1, 0]


def func(x, k, m):
    return k * x + m


def find_foto_tvärsnitt(energi_list, foto_list, energi):
    energi_list = np.array(energi_list)
    foto_list = np.array(foto_list)

    diff = np.abs(energi_list - energi)
    closest_indices = np.argsort(diff)[:3]
    foto_close = foto_list[closest_indices]
    energi_close = energi_list[closest_indices]
    # print('energi close: ', energi_close)
    print('foto close: ', foto_close)

    # linjär anpassning
    popt, _ = curve_fit(func, energi_close, foto_close)

    k, m = popt
    print('k, m: ',k, m)
    foto_target = func(energi, k, m)
    print(foto_target)
    return foto_target

energi = 10**5
foto_target = find_foto_tvärsnitt(energi_list, foto_list, energi)

x_data = [energi_list, energi]
y_data = [foto_list, foto_target]
scatter = [2, 1]
label_data = 'foto tvärsnitt'
marker = ['o', 'X']
color = ['blue', 'red']

fig = plot_stuff(x_data, y_data, scatter, label_data,
                 marker, color, x_label='energi (eV)', y_label='tvärsnitt (cm^2)', title='1',
                 fig_size=(10, 10), symbol_size=100, font_size=30, alpha=1, line_width=2, x_lim=(0, 0), y_lim=(0, 0),
                 grid=False, x_scale='log', y_scale='log')

fig.savefig('foto tvärsnitt', bbox_inches='tight')

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

# def växelverkan(energi):
#
