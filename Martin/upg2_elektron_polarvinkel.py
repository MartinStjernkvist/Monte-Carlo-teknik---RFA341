from imports import *


def elektron_theta_ny(elektron_energi, scatterpower_data, rho_medium):  # Energi i elektron i tabell är i eV

    # Anta att tumören är i en vävnad alltså se på vatten för tvärsnitten

    """
    #Läser in excel filen och tar ut koloumn med energi och scatter power
    file_scatter=pd.read_excel(r'MC_Linnea/Monte Carlo, test.xlsx', sheet_name='Scattering power') #Behövde skriva av en egen Excel fil från ICRU 35
    Energi_data=file_scatter['E/MeV'].to_list()
    Mass_Scatterpower_data=file_scatter['Scatter power Water '].to_list()  #T/rho för vatten
    """

    # Omvandlar från MeV till eV och till en np.array
    energi_MeV_list = scatterpower_data[:, 0]
    energi_list = np.array([(lambda x: x * 10 ** 6)(x) for x in energi_MeV_list])

    mass_scatterpower_list = np.array(scatterpower_data[:, 1])

    # Hittar närmast energi som liknar elektron energin och ta fram scatter power:

    # Tar index för närmaste energi på elektronen
    diff = np.abs(energi_list - elektron_energi)
    closest_indices = np.argsort(diff)[:2]

    # Få scattering power nära energierna och energivärdena närmast
    T_close = mass_scatterpower_list[closest_indices] * rho_medium * 10 ** (-3)  # kg/m^3 till g/cm^3
    energi_close = energi_list[closest_indices]

    # Linjär interpolera och få fram theta_s
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        T = T_close[0]

    else:
        T = T_close[0] + (elektron_energi - energi_close[0]) * (
                T_close[1] - T_close[0]) / (
                           energi_close[1] - energi_close[0])

    R = np.random.random()  # Slumpmässig tal mellan 0-1

    dl = 10**(-20)

    theta_ny = np.sqrt(-T * dl * np.log(1 - R))

    theta_ny = random.choice([1, -1]) * theta_ny
    print('theta ny', theta_ny)

    return theta_ny
