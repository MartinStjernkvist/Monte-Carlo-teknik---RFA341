from imports import *


def scattering_power_från_energi(elektron_energi, scatterpower_data, rho_medium):
    """
    Funktion som linjärinterpolerar scattering power utifrån energi.
    :param elektron_energi: Elektronens energi.
    :param scatterpower_data: Scattering power tabelldata.
    :param rho_medium: Mediumets densitet (kg/m^3).
    :return: Scattering power (rad^2 / cm).
    """

    # Omvandlar från MeV till eV och till en np.array
    energi_MeV_list = scatterpower_data[:, 0]
    energi_list = np.array([(lambda x: x * 10 ** 6)(x) for x in energi_MeV_list])

    mass_scattering_power_list = np.array(scatterpower_data[:, 1])

    # Hittar närmast energi som liknar elektron energin och ta fram scatter power:

    # Tar index för närmaste energi på elektronen.
    diff = np.abs(energi_list - elektron_energi)
    closest_indices = np.argsort(diff)[:2]

    # Få scattering power nära energierna och energivärdena närmast.
    T_close_cm = mass_scattering_power_list[closest_indices] * rho_medium * 10 ** (-3)  # Densitet: kg/m^3 till g/cm^3
    T_close = T_close_cm * 10 ** 2  # scattering power / m

    energi_close = energi_list[closest_indices]

    # Linjärinterpolering av scattering power.
    if energi_close[1] - energi_close[0] < 10 ** (-15):
        T = T_close[0]

    else:
        T = T_close[0] + (elektron_energi - energi_close[0]) * (
                T_close[1] - T_close[0]) / (
                    energi_close[1] - energi_close[0])

    return T
