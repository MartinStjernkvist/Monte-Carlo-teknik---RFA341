from imports import *

from sampla_energi_start import energi_start, Lu177_energi, Lu177_intensitet, Lu177_sannolikhet
from matriser import slicad_fantom_matris, slicad_njure_matris, slicad_benmärg_matris
from sampla_position_start import position_start
from sampla_riktning_och_steg_start import riktning_start, första_steg
from sampla_steglängd import medelvägslängd
from sampla_växelverkan import växelverkan, växelverkan_slimmad
from transformation_3d import transformera_koordinatsystem
from attenueringsdata import attenueringsdata


def run_MC(iterationer, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris, slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi, radionuklid_intensitet, radionuklid_sannolikhet):

    df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)
    df_anatomidefinitioner = pd.read_excel(anatomidefinitioner_file, index_col=None)
    df_tvärsnitt = pd.read_excel(tvärsnitt_file, index_col=None)

    x_size, y_size, z_size = slicad_fantom_matris.shape
    benmärg_matris_deponerad_energi = np.zeros((x_size, y_size, z_size))

    utanför_fantom = 0
    vxv_foto = 0
    # vxv_compton = 0
    # vxv_rayleigh = 0
    träff = 0

    for i in range(iterationer):
        attenuerad = 0
        # print(i)

        # start: sampla position, riktning och energi
        foton_energi = energi_start(radionuklid_energi, radionuklid_intensitet, radionuklid_sannolikhet)
        x_start, y_start, z_start = position_start(slicad_njure_matris)
        theta, phi = riktning_start()

        voxel_värde = slicad_fantom_matris[x_start, y_start, z_start]
        instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
        mu = instans.mu()

        # steglängd: sampla medelvägslängden från inverstransformerad attenueringsfunktion
        steglängd = medelvägslängd(mu)
        # print(f'steglängd: {steglängd}')

        # steg: gå steget till ny position i startriktning
        x, y, z = första_steg(theta, phi, steglängd, x_start, y_start, z_start)
        x_round, y_round, z_round = round(x), round(y), round(z)

        if (
                x_round < 0
                or x_round >= x_size
                or y_round < 0
                or y_round >= y_size
                or z_round < 0
                or z_round >= z_size
        ):
            utanför_fantom += 1
            i += 1

        else:
            while attenuerad == 0:

                voxel_värde = slicad_fantom_matris[x_round, y_round, z_round]
                instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
                mu = instans.mu()

                if slicad_fantom_matris[x_round, y_round, z_round] == 0:
                    utanför_fantom += 1
                    attenuerad = 1

                else:
                    instans = växelverkan_slimmad(foton_energi, df_tvärsnitt)
                    vxv = instans.bestäm_växelverkan_slimmad()
                    # print(f'energi: {foton_energi * 10 ** (-3)} keV, vxv: {vxv}')

                    if vxv == 'foto':
                        vxv_foto += 1

                        if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                            träff += 1

                            benmärg_matris_deponerad_energi[x_round, y_round, z_round] += foton_energi

                            print(
                                f'foto: {energideponering_compton * 10 ** (-3):.2f} keV i voxel [{round(x), round(y), round(z)}]')

                        attenuerad = 1

                    elif vxv == 'compton':
                        # vxv_compton += 1

                        # sampla energideponering, vinkel
                        # theta_compton = sampla_theta_compton
                        # mu_nytt energiberoende?
                        theta_compton = 1
                        phi_compton = 2 * pi * random.rand()

                        energideponering_compton = foton_energi * 0.5  # placeholder

                        if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                            träff += 1

                            benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energideponering_compton

                            print(
                                f'compton: {energideponering_compton * 10 ** (-3):.2f} keV i voxel [{round(x), round(y), round(z)}]')

                        foton_energi = foton_energi - energideponering_compton

                        steglängd_compton = medelvägslängd(mu)
                        vektor_compton, dx_compton, dy_compton, dz_compton = transformera_koordinatsystem(steglängd,
                                                                                                          phi,
                                                                                                          theta,
                                                                                                          steglängd_compton,
                                                                                                          phi_compton,
                                                                                                          theta_compton)

                        x_round = round(x + dx_compton / voxel_sidlängd)
                        y_round = round(y + dy_compton / voxel_sidlängd)
                        z_round = round(z + dz_compton / voxel_sidlängd)

                        if (
                                x_round < 0
                                or x_round >= x_size
                                or y_round < 0
                                or y_round >= y_size
                                or z_round < 0
                                or z_round >= z_size
                        ):
                            utanför_fantom += 1
                            attenuerad = 1

                        theta, phi = theta_compton, phi_compton
                        steglängd = steglängd_compton

                    elif vxv == 'rayleigh':
                        # vxv_rayleigh += 1

                        theta_rayleigh = 1
                        phi_rayleigh = 2 * pi * random.rand()
                        steglängd_rayleigh = medelvägslängd(mu)

                        vektor_rayleigh, dx_rayleigh, dy_rayleigh, dz_rayleigh = transformera_koordinatsystem(
                            steglängd,
                            phi,
                            theta,
                            steglängd_rayleigh,
                            phi_rayleigh,
                            theta_rayleigh)

                        x_round = round(x + dx_rayleigh / voxel_sidlängd)
                        y_round = round(y + dy_rayleigh / voxel_sidlängd)
                        z_round = round(z + dz_rayleigh / voxel_sidlängd)

                        if (
                                x_round < 0
                                or x_round >= x_size
                                or y_round < 0
                                or y_round >= y_size
                                or z_round < 0
                                or z_round >= z_size
                        ):
                            utanför_fantom += 1
                            attenuerad = 1

                        theta, phi = theta_rayleigh, phi_rayleigh
                        steglängd = steglängd_rayleigh

    print(f'max värdet av matrisen: {np.max(benmärg_matris_deponerad_energi)}')
    print(f'utanför: {utanför_fantom}')
    print(f'foto: {vxv_foto}')
    # print(f'rayleigh: {vxv_rayleigh}')
    # print(f'compton: {vxv_compton}')
    print(f'träffar: {träff}')
    return benmärg_matris_deponerad_energi


#   ----------------------------------------------------------------------
#   KÖR KODEN
#   ----------------------------------------------------------------------
if __name__ == "__main__":
    start = time.time()

    iterationer = 10 **2

    radionuklid_energi = Lu177_energi
    radionuklid_intensitet = Lu177_intensitet
    radionuklid_sannolikhet = Lu177_sannolikhet

    benmärg_matris_deponerad_energi = run_MC(iterationer, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris, slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi, radionuklid_intensitet, radionuklid_sannolikhet)

    end_time(start)

    np.save('benmärg_matris_deponerad_energi.npy', benmärg_matris_deponerad_energi)
