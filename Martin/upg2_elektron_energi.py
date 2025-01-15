from imports import *

#   -----------------------------------
#   Sönderfallsdata för Y-90.
#   -----------------------------------

file_Y90 = pd.read_excel(Y90_file)

# Ta ut mätpunkter med intensitet och energi.
Energi_Y90 = file_Y90['Energy (MeV)']  # MeV
Intensitet_Y90 = file_Y90['#/nt']


#   -----------------------------------
#   Kurvanpassning till sönderfallsdatan.
#   -----------------------------------

def polynom_funktion(x, a, b, c, d, e, f):
    """
    5:gradigt polynom för kurvanpassning.
    """
    return a * x ** 5 + b * x ** 4 + c * x ** 3 + d * x ** 2 + e * x + f

# Interpolera.
params, cv = curve_fit(polynom_funktion, Energi_Y90, Intensitet_Y90)
a, b, c, d, e, f = params

olika_energier = np.linspace(np.min(Energi_Y90), np.max(Energi_Y90), 10_000)


#   -----------------------------------
#   Använd rejektionsmetoden för att sampla elektronenergi.
#   -----------------------------------
def elektron_energi_start():
    """
    Funktion som samplar elektronens startenergi.
    Använder rejektionsmetoden.
    """

    # Initiera värden.
    elektron_energi = 0
    hittat = 0

    # Medan inget värde accepterats:
    while hittat == 0:

        # Slumpmässiga värden på energi och intensitet.
        x_rand = np.random.rand() * np.max(Energi_Y90)
        y_rand = np.random.rand() * np.max(Intensitet_Y90)

        # Rejektionsmetoden för att hitta en energi.
        if y_rand < polynom_funktion(x_rand, *params):

            # Om hittat: erhåll energin och avsluta while loopen.
            hittat = 1
            elektron_energi = x_rand

        else:
            hittat = 0

    return elektron_energi


if __name__ == "__main__":
    plt.plot(olika_energier, polynom_funktion(olika_energier, *params))

    plt.scatter(Energi_Y90, Intensitet_Y90)

    # Visa figuren
    plt.show()
