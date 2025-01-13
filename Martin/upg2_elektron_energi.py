from imports import *

# Slumpa ut energi på betakällan


# Hittar filen och ta info från excelfilen
# Läser en Excel fil

file_Y90 = pd.read_excel(Y90_file)

# Ta ut värderna på energin och intensiteten betakällan
Energi_Y90 = file_Y90['Energy (MeV)']  # MeV
Intensitet_Y90 = file_Y90['#/nt']
# Plottar ut värdena från excel filen

plt.scatter(Energi_Y90, Intensitet_Y90)


# Plottar ut punkterna i excelfilen och gör en kurvanpassning

# def polynom_funktion(x, a, b, c, d):
#     return a * x ** 3 + b * x ** 2 + c * x + d

def polynom_funktion(x, a, b, c, d, e, f):
    return a * x ** 5 + b * x ** 4 + c * x ** 3 + d * x ** 2 + e * x + f


params, cv = curve_fit(polynom_funktion, Energi_Y90, Intensitet_Y90)
# print(*params)
# a, b, c, d = params
a, b, c, d, e, f = params
olika_energier = np.linspace(np.min(Energi_Y90), np.max(Energi_Y90), 1000)

plt.plot(olika_energier, polynom_funktion(olika_energier, *params))

# Visa figuren
plt.show()

f_max = np.max(polynom_funktion(olika_energier, *params))


# print(f_max) #Test om det stämmer


# Använd Rejektionsmetoden för att sampla elektronenergi


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
        x_sampel = np.random.random() * Skärpunkt_0
        if np.random.random() <= polynom_funktion(x_sampel, *params) / f_max:
            Elektron_energi = x_sampel
            return Elektron_energi
        else:
            continue
