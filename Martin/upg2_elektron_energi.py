from imports import *

#   ----------------------------------------------------------------------
#   Sönderfallsdata för Y-90.
#   ----------------------------------------------------------------------

file_Y90 = pd.read_excel(Y90_file)

# Ta ut mätpunkter med intensitet och energi.
Energi_Y90 = file_Y90['Energy (MeV)']  # MeV
Intensitet_Y90 = file_Y90['#/nt']


#   ----------------------------------------------------------------------
#   Kurvanpassning till sönderfallsdatan.
#   ----------------------------------------------------------------------

def polynom_funktion(x, a, b, c, d, e, f):
    return a * x ** 5 + b * x ** 4 + c * x ** 3 + d * x ** 2 + e * x + f


params, cv = curve_fit(polynom_funktion, Energi_Y90, Intensitet_Y90)
a, b, c, d, e, f = params

olika_energier = np.linspace(np.min(Energi_Y90), np.max(Energi_Y90), 10_000)

f_max = np.max(polynom_funktion(olika_energier, *params))


#   ----------------------------------------------------------------------
#   Använd rejektionsmetoden för att sampla elektronenergi.
#   ----------------------------------------------------------------------

# Tar fram värdet närmast skärningspunkten i x-axeln
def närmast(lista, tal):
    lista = np.array(lista)
    närmsta_index = (np.abs(lista - tal)).argmin()
    return lista[närmsta_index], närmsta_index


close, _ = närmast(polynom_funktion(olika_energier, *params), 0)

for i in range(len(olika_energier)):
    if polynom_funktion(olika_energier[i], *params) == close:
        Skärpunkt_0 = olika_energier[i]
    else:
        continue


# Sampla ett x värde mellan 0 och skärningspunkten för att få elektronens energi

def elektron_startenergi():
    while True:
        x_sampel = np.random.rand() * Skärpunkt_0
        if np.random.rand() <= polynom_funktion(x_sampel, *params) / f_max:
            Elektron_energi = x_sampel
            return Elektron_energi
        else:
            continue


def elektron_energi_start():

    hittat = 0
    while hittat == 0:
        x_rand = np.random.rand() * np.max(Energi_Y90)
        y_rand = np.random.rand() * np.max(Intensitet_Y90)

        if y_rand < polynom_funktion(x_rand, *params):
            hittat = 1

        else:
            hittat = 0

    return elektron_energi_start


if __name__ == "__main__":
    plt.plot(olika_energier, polynom_funktion(olika_energier, *params))

    plt.scatter(Energi_Y90, Intensitet_Y90)

    # Visa figuren
    plt.show()
