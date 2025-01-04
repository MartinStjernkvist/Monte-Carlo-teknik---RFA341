from imports import *

from sampla_energi_start import energi_start, Lu177_energi, Lu177_intensitet, Lu177_sannolikhet
from matriser import slicad_fantom_matris, slicad_njure_matris, slicad_benmärg_matris
from sampla_position_start import position_start
from sampla_riktning_och_steg_start import riktning_koherent, steg
from sampla_steglängd import medelvägslängd
from sampla_växelverkan import växelverkan
from transformation_3d import ny_transformera_koordinatsystem
from attenueringsdata import attenueringsdata
from sampla_compton_vinkel import compton_vinkel_och_energiförlust
from sampla_foto_vxv import foto_vxv
from bestäm_om_attenuerad import bestäm_om_attenuerad
from input_upg1_multiprocess import iterationer_tot, antal_cores, iterationer_dummy


def run_MC_multiprocess(args):
    (start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
     slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi, radionuklid_intensitet,
     radionuklid_sannolikhet) = args

    #   ----------------------------------------------------------------------
    #   Läs in data
    #   ----------------------------------------------------------------------
    df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)
    df_anatomidefinitioner = pd.read_excel(anatomidefinitioner_file, index_col=None)
    df_tvärsnitt = pd.read_excel(tvärsnitt_file, index_col=None)

    # Matrisdimensionerna, viktigt för att bedöma ifall fotonen är innanför matrisen eller inte.
    x_size, y_size, z_size = slicad_fantom_matris.shape
    benmärg_matris_deponerad_energi = np.zeros((x_size, y_size, z_size))

    #   ----------------------------------------------------------------------
    #   Håll reda på vad som händer när koden körs
    #   ----------------------------------------------------------------------
    utanför_fantom = 0
    vxv_foto = 0
    träff_benmärg = 0

    #   ----------------------------------------------------------------------
    #   Loopar igenom alla iterationer
    #   ----------------------------------------------------------------------
    for i in range(start, end):

        # initiera loopen med att attenuerad = 0
        attenuerad = 0
        # print(i)

        # start: sampla position, riktning och energi
        foton_energi = energi_start(radionuklid_energi, radionuklid_intensitet, radionuklid_sannolikhet)
        x_start, y_start, z_start = position_start(slicad_njure_matris)
        theta, phi = riktning_koherent()

        voxel_värde = slicad_fantom_matris[x_start, y_start, z_start]
        instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
        mu = instans.mu()

        # Sampla medelvägslängden från inverstransformerad attenueringsfunktion.
        steglängd = medelvägslängd(mu)
        # print(f'steglängd: {steglängd}')

        # Gå steget till ny position från startpositionen i startriktningen.
        x, y, z = steg(theta, phi, steglängd, x_start, y_start, z_start)
        x_round, y_round, z_round = round(x), round(y), round(z)

        #   ----------------------------------------------------------------------
        #   Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen
        #   Ekvationen förekommer nedan efter varje nytt steg tas, dock utan if-statement.
        #   ----------------------------------------------------------------------
        attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size, utanför_fantom, slicad_fantom_matris, foton_energi)
        if attenuerad == 1:
            i += 1

        else:

            #   ----------------------------------------------------------------------
            #   Loopa under tiden som fotonen inte attenuerats och fortfarande är i matrisen.
            #   ----------------------------------------------------------------------
            while attenuerad == 0:

                # Identifiera vilken voxel fotonen befinner sig i.
                voxel_värde = slicad_fantom_matris[x_round, y_round, z_round]
                instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
                mu = instans.mu()

                # Bestäm vilken typ av växelverkan som sker vid nya positionen.
                instans = växelverkan(foton_energi, df_tvärsnitt)
                vxv = instans.bestäm_växelverkan()
                # print(f'energi: {foton_energi * 10 ** (-3)} keV, vxv: {vxv}')

                #   ----------------------------------------------------------------------
                #
                #   Fotoabsorption.
                #
                #   ----------------------------------------------------------------------
                if vxv == 'foto':
                    vxv_foto += 1

                    # Bestäm om fluorescens sker eller inte.
                    energi_deponering, attenuerad = foto_vxv(foton_energi)
                    foton_energi = foton_energi - energi_deponering

                    # Registrera endast energideponeringen om nuvarande voxel är i benmärg matrisen.
                    if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                        träff_benmärg += 1

                        benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energi_deponering

                        print(
                            f'foto: {energi_deponering:.2f} eV i voxel [{round(x), round(y), round(z)}]')

                    # Om fluorescens sker -> följ ny foton.
                    if attenuerad == 0:

                        # Sampla spridningsvinklar
                        theta_foto, phi_foto = riktning_koherent()

                        voxel_värde = slicad_fantom_matris[x_round, y_round, z_round]
                        instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata,
                                                   df_anatomidefinitioner)
                        mu = instans.mu()
                        steglängd_foto = medelvägslängd(mu)

                        x, y, z = steg(theta, phi, steglängd, x_round, y_round, z_round)
                        x_round, y_round, z_round = round(x), round(y), round(z)

                        # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                        attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size, utanför_fantom, slicad_fantom_matris, foton_energi)

                        # Ingångsvärden till koordinat-transformeringen som behöver genomföras ifall växelverkan = compton eller rayleigh.
                        theta, phi = theta_foto, phi_foto
                        steglängd = steglängd_foto

                #   ----------------------------------------------------------------------
                #
                #   Comptonspridning.
                #
                #   ----------------------------------------------------------------------
                elif vxv == 'compton':

                    # Sampla spridningsvinklar, samt energideponeringen vid Comptonväxelverkan.
                    theta_compton, foton_energi, energideponering_compton = compton_vinkel_och_energiförlust(
                        foton_energi)
                    phi_compton = 2 * pi * random.rand()

                    if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                        träff_benmärg += 1

                        benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energideponering_compton

                        print(
                            f'compton: {energideponering_compton:.2f} eV i voxel [{round(x), round(y), round(z)}]')

                    steglängd_compton = medelvägslängd(mu)

                    # Koordinattransformation, eftersom spridningsvinklarna för Comptonspridning inte kan samplas uniformt.
                    dx_compton, dy_compton, dz_compton = ny_transformera_koordinatsystem(steglängd, phi, theta, steglängd_compton, phi_compton, theta_compton)

                    x = x + dx_compton
                    y = y + dy_compton
                    z = z + dz_compton

                    x_round, y_round, z_round = round(x), round(y), round(z)

                    # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                    attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size, utanför_fantom, slicad_fantom_matris, foton_energi)

                    # Ingångsvärden till koordinat-transformeringen (om nästa växelverkan är Comptonspridning eller Rayleighspridning).
                    theta, phi = theta_compton, phi_compton
                    steglängd = steglängd_compton

                #   ----------------------------------------------------------------------
                #
                #   Rayleighspridning
                #
                #   ----------------------------------------------------------------------
                elif vxv == 'rayleigh':
                    # vxv_rayleigh += 1

                    theta_rayleigh = 1 # ERSÄTT MED THOMPSON TVÄRSNITT

                    voxel_värde = slicad_fantom_matris[x_round, y_round, z_round]
                    instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata,
                                               df_anatomidefinitioner)
                    mu = instans.mu()

                    phi_rayleigh = 2 * pi * np.random.rand()
                    steglängd_rayleigh = medelvägslängd(mu)

                    # Koordinattransformation, eftersom spridningsvinklarna för Rayleighspridning inte kan samplas uniformt.
                    dx_rayleigh, dy_rayleigh, dz_rayleigh = ny_transformera_koordinatsystem(
                        steglängd,
                        phi,
                        theta,
                        steglängd_rayleigh,
                        phi_rayleigh,
                        theta_rayleigh)

                    x = x + dx_rayleigh
                    y = y + dy_rayleigh
                    z = z + dz_rayleigh

                    x_round, y_round, z_round = round(x), round(y), round(z)

                    # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                    attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size, utanför_fantom, slicad_fantom_matris, foton_energi)

                    # Ingångsvärden till koordinat-transformeringen (om nästa växelverkan är Comptonspridning eller Rayleighspridning).
                    theta, phi = theta_rayleigh, phi_rayleigh
                    steglängd = steglängd_rayleigh

    print(f'max värdet av matrisen: {np.max(benmärg_matris_deponerad_energi)}')
    print(f'utanför: {utanför_fantom}')
    print(f'foto: {vxv_foto}')
    # print(f'rayleigh: {vxv_rayleigh}')
    # print(f'compton: {vxv_compton}')
    print(f'träffar: {träff_benmärg}')
    return benmärg_matris_deponerad_energi


# def input_viktiga_parametrar():
#     root = tk.Tk()
#     root.withdraw()
#
#     tot_iterationer = 10 ** (simpledialog.askfloat("Iterationer", "Ange magnitud:"))
#     antal_cores = round(int(simpledialog.askfloat("Cores", "[VIKTIGT] Ange antal kärnor:")))
#
#     if antal_cores >= 8:
#         antal_cores = 1
#         print('FÖR STORT VÄRDE PÅ ANTAL CORES')
#
#     if tot_iterationer >= 10**7:
#         tot_iterationer = 1
#         print('FÖR STORT VÄRDE PÅ ANTAL ITERATIONER')
#
#     return tot_iterationer, antal_cores

#   ----------------------------------------------------------------------
#   KÖR KODEN
#   ----------------------------------------------------------------------
if __name__ == "__main__":

    radionuklid_energi = Lu177_energi
    radionuklid_intensitet = Lu177_intensitet
    radionuklid_sannolikhet = Lu177_sannolikhet

    #   ----------------------------------------------------------------------
    #   DUMMY RUN
    #   ----------------------------------------------------------------------
    print(
        '\n----------------------------------------------------------------------\nDUMMY RUN\n----------------------------------------------------------------------\n')
    start = time.time()

    chunk_storlek = iterationer_dummy // antal_cores
    chunk_ranges = [(i * chunk_storlek, (i + 1) * chunk_storlek) for i in range(antal_cores)]
    chunk_ranges[-1] = (chunk_ranges[-1][0], iterationer_dummy)

    args_packed = [(start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
                    slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi,
                    radionuklid_intensitet, radionuklid_sannolikhet) for start, end in chunk_ranges]

    with mp.Pool(antal_cores) as pool:
        partial_results = pool.map(run_MC_multiprocess, args_packed)

    benmärg_matris_deponerad_energi = np.sum(partial_results, axis=0)

    end_time(start)

    #   ----------------------------------------------------------------------
    #   ACTUAL
    #   ----------------------------------------------------------------------

    print(
        '\n----------------------------------------------------------------------\nACTUAL RUN\n----------------------------------------------------------------------\n')
    start = time.time()

    radionuklid_energi = Lu177_energi
    radionuklid_intensitet = Lu177_intensitet
    radionuklid_sannolikhet = Lu177_sannolikhet

    chunk_storlek = iterationer_tot // antal_cores
    chunk_ranges = [(i * chunk_storlek, (i + 1) * chunk_storlek) for i in range(antal_cores)]
    chunk_ranges[-1] = (chunk_ranges[-1][0], iterationer_tot)

    args_packed = [(start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
                    slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi,
                    radionuklid_intensitet, radionuklid_sannolikhet) for start, end in chunk_ranges]

    with mp.Pool(antal_cores) as pool:
        partial_results = pool.map(run_MC_multiprocess, args_packed)

    benmärg_matris_deponerad_energi = np.sum(partial_results, axis=0)

    end_time(start)

    np.save('resultat_multiprocess.npy', benmärg_matris_deponerad_energi)

# 24.4 seconds
