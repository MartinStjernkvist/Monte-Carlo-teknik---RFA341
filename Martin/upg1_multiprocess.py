from imports import *

from upg1_sampla_energi_start import energi_start
from upg1_matriser import slicad_fantom_matris, slicad_njure_matris, slicad_benmärg_matris
from upg1_sampla_position_start import position_start
from upg1_sampla_riktning_och_steg_start import riktning_uniform, steg
from upg1_sampla_steglängd import medelvägslängd
from upg1_sampla_växelverkan import växelverkan
from upg1_transformation_3d import steg_transformera_koordinatsystem_3d
from upg1_attenueringsdata import attenueringsdata
from upg1_sampla_compton import compton_vinkel_och_energiförlust
from upg1_sampla_foto_vxv import foto_vxv
from upg1_bestäm_om_attenuerad import bestäm_om_attenuerad
from upg12_steg_transformation import ny_steg_transformera_koordinatsystem_3d
from upg12_förflyttning import förflyttning


# from input_upg1_multiprocess import iterationer_tot, antal_cores, iterationer_dummy


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

    # Skapa tom matris, för att sedan fylla på med energideponering i voxklarna.
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

        # Initiera loopen med att attenuerad = 0.
        attenuerad = 0
        # print(i)

        # Start: sampla position, riktning och energi
        foton_energi = energi_start(radionuklid_energi, radionuklid_sannolikhet)
        x_start, y_start, z_start = position_start(slicad_njure_matris)
        theta, phi = riktning_uniform()

        voxel_värde = slicad_fantom_matris[x_start, y_start, z_start]
        instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
        mu = instans.mu_max()

        # Sampla medelvägslängden från inverstransformerad attenueringsfunktion.
        steglängd = medelvägslängd(mu)
        # print(f'steglängd: {steglängd}')

        # Gå steget till ny position från startpositionen i startriktningen.
        # Eftersom voxlarna har diskreta positioner måste avrundning till närmaste heltal göras.
        dx, dy, dz = steg(theta, phi, steglängd)
        x, y, z, x_round, y_round, z_round = förflyttning(x_start, y_start, z_start, dx, dy, dz)

        #   ----------------------------------------------------------------------
        #   Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen
        #   Ekvationen förekommer nedan efter varje nytt steg tas, dock utan if-statement.
        #   ----------------------------------------------------------------------

        attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size,
                                                          utanför_fantom, slicad_fantom_matris, foton_energi)
        # Kanske kolla om den fortsätter eller vxv?

        if attenuerad == 1:
            i += 1

        else:

            #   ----------------------------------------------------------------------
            #   Loopa under tiden som fotonen inte attenuerats och fortfarande är i matrisen.
            #   ----------------------------------------------------------------------
            while attenuerad == 0:
                # print('steglängd: ', steglängd)

                # Identifiera vilken voxel fotonen befinner sig i.
                voxel_värde = slicad_fantom_matris[x_round, y_round, z_round]
                instans = attenueringsdata(voxel_värde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
                mu = instans.mu_max()

                # Bestäm vilken typ av växelverkan som sker vid nya positionen.
                instans = växelverkan(foton_energi, df_tvärsnitt)
                vxv = instans.bestäm_växelverkan()
                # print(f'energi: {foton_energi * 10 ** (-3)} keV, vxv: {vxv}')

                #   ----------------------------------------------------------------------
                #   Fotoabsorption.
                #   ----------------------------------------------------------------------
                if vxv == 'foto':
                    vxv_foto += 1

                    # Bestäm om fluorescens sker eller inte.
                    energi_deponering, attenuerad = foto_vxv(foton_energi)
                    foton_energi = foton_energi - energi_deponering

                    # Registrera energideponering ifall fotonen växelverkar i en voxel med benmärg.
                    if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                        träff_benmärg += 1

                        benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energi_deponering

                        print(
                            f'foto: {energi_deponering:.0f} eV i voxel [{round(x), round(y), round(z)}]')

                    # Om fluorescens sker -> följ ny foton.
                    if attenuerad == 0:
                        # Sampla spridningsvinklar (uniformt samplade).
                        theta_foto, phi_foto = riktning_uniform()

                        # Ta ett nytt steg.
                        steglängd_foto = medelvägslängd(mu)
                        dx, dy, dz = steg(theta, phi, steglängd)
                        x, y, z, x_round, y_round, z_round = förflyttning(x_start, y_start, z_start, dx, dy, dz)

                        # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                        attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size,
                                                                          z_size, utanför_fantom, slicad_fantom_matris,
                                                                          foton_energi)

                        # Ingångsvärden till koordinat-transformeringen som behöver genomföras ifall växelverkan = compton eller rayleigh.
                        theta, phi = theta_foto, phi_foto
                        steglängd = steglängd_foto

                #   ----------------------------------------------------------------------
                #   Comptonspridning.
                #   ----------------------------------------------------------------------
                elif vxv == 'compton':

                    # Sampla spridningsvinklar, samt energideponeringen vid Comptonväxelverkan.
                    theta_compton, foton_energi, energideponering_compton = compton_vinkel_och_energiförlust(
                        foton_energi)
                    phi_compton = 2 * pi * np.random.rand()

                    # Registrera energideponering ifall fotonen växelverkar i en voxel med benmärg.
                    if slicad_benmärg_matris[x_round, y_round, z_round] != 0:
                        träff_benmärg += 1

                        benmärg_matris_deponerad_energi[x_round, y_round, z_round] += energideponering_compton

                        print(
                            f'compton: {energideponering_compton:.0f} eV i voxel [{round(x), round(y), round(z)}]')

                    steglängd_compton = medelvägslängd(mu)

                    # Koordinattransformation, eftersom spridningsvinklarna för Comptonspridning inte kan samplas uniformt.
                    dx_compton, dy_compton, dz_compton = ny_steg_transformera_koordinatsystem_3d(steglängd, phi, theta,
                                                                                                 steglängd_compton,
                                                                                                 phi_compton,
                                                                                                 theta_compton)

                    # Ta ett nytt steg.
                    x, y, z, x_round, y_round, z_round = förflyttning(x, y, z, dx_compton, dy_compton, dz_compton)

                    # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                    attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size,
                                                                      utanför_fantom, slicad_fantom_matris,
                                                                      foton_energi)

                    # Ingångsvärden till koordinat-transformeringen (om nästa växelverkan är Comptonspridning eller Rayleighspridning).
                    theta, phi = theta_compton, phi_compton
                    steglängd = steglängd_compton

                #   ----------------------------------------------------------------------
                #   Rayleighspridning.
                #   ----------------------------------------------------------------------
                elif vxv == 'rayleigh':

                    # Sampla spridningsvinklar och steglängd.
                    theta_rayleigh = np.arccos(-1 + 2 * np.random.rand())  # ERSÄTT MED THOMSON TVÄRSNITT
                    phi_rayleigh = 2 * pi * np.random.rand()
                    steglängd_rayleigh = medelvägslängd(mu)

                    # Koordinattransformation, eftersom spridningsvinklarna för Rayleighspridning inte kan samplas uniformt.
                    dx_rayleigh, dy_rayleigh, dz_rayleigh = ny_steg_transformera_koordinatsystem_3d(
                        steglängd,
                        phi,
                        theta,
                        steglängd_rayleigh,
                        phi_rayleigh,
                        theta_rayleigh)

                    # Ta ett nytt steg.
                    x, y, z, x_round, y_round, z_round = förflyttning(x, y, z, dx_rayleigh, dy_rayleigh, dz_rayleigh)

                    # Om foton hamnar utanför fantommatrisen -> kasta ut foton ur loopen.
                    attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size,
                                                                      utanför_fantom, slicad_fantom_matris,
                                                                      foton_energi)

                    # Ingångsvärden till koordinat-transformeringen (om nästa växelverkan är Comptonspridning eller Rayleighspridning).
                    theta, phi = theta_rayleigh, phi_rayleigh
                    steglängd = steglängd_rayleigh

    # Några print statements, för att hålla reda på vad som händer när koden kör.
    print(f'max värdet av matrisen: {np.max(benmärg_matris_deponerad_energi)}')
    print(f'utanför: {utanför_fantom}')
    print(f'foto: {vxv_foto}')
    print(f'träffar: {träff_benmärg}')

    return benmärg_matris_deponerad_energi


if __name__ == "__main__":
    radionuklid_energi = Lu177_energi
    radionuklid_intensitet = Lu177_intensitet
    radionuklid_sannolikhet = Lu177_sannolikhet


    def inputs_riktig_körning():
        print('Standard eller inte?')
        input_standard = input('Om standard: s, Annars: vad som helst: ')

        if input_standard == 's':
            antal_cores = 8
            iterationer_dummy = 10 ** 3
            iterationer_tot = 10 ** 5

        else:
            print(
                '\n----------------------------------------------------------------------\nVIKTIGT:\n----------------------------------------------------------------------\nAnge antal processor kärnor')
            input_antal_cores = input('Antal kärnor: ')

            if eval(input_antal_cores) > 8:
                antal_cores = 1
            else:
                antal_cores = eval(input_antal_cores)

            print(
                '\n----------------------------------------------------------------------\nDUMMY:\n----------------------------------------------------------------------\nAnge magnitud: ex 3 -> 10^3 iterationer')
            input_dummy_magnitud_iterationer = input('Magnitud: ')

            if eval(input_dummy_magnitud_iterationer) > 4:
                iterationer_dummy = 10 ** 3
            else:
                iterationer_dummy = 10 ** (eval(input_dummy_magnitud_iterationer))

            print(
                '\n----------------------------------------------------------------------\nRIKTIG:\n----------------------------------------------------------------------\nAnge skalär och magnitud: ex 5 och 5 -> 5 * 10^5 iterationer')
            input_riktig_skalär_iterationer = input('Skalär: ')
            input_riktig_magnitud_iterationer = input('Magnitud: ')

            if eval(input_riktig_magnitud_iterationer) >= 8:
                iterationer_tot = 10 ** 3
            else:
                iterationer_tot = eval(input_riktig_skalär_iterationer) * 10 ** (
                    eval(input_riktig_magnitud_iterationer))

        print('antal_cores, iterationer_dummy, iterationer_tot: ', antal_cores, iterationer_dummy, iterationer_tot)
        return antal_cores, iterationer_dummy, iterationer_tot


    antal_cores, iterationer_dummy, iterationer_tot = inputs_riktig_körning()

    dictionary = {
        "antal_cores": antal_cores,
        "iterationer_dummy": iterationer_dummy,
        "iterationer_tot": iterationer_tot
    }

    json_object = json.dumps(dictionary)

    with open('inputs_upg1_multiprocess.json', 'w') as f:
        f.write(json_object)
        f.close()

    #   ----------------------------------------------------------------------
    #   Dummy run - för att snabba på den riktiga körningen av koden.
    #   ----------------------------------------------------------------------
    print(
        '\n----------------------------------------------------------------------\nDUMMY RUN\n----------------------------------------------------------------------\n')
    start = time.time()

    # Lite kod för att dela upp arbetet i flera processer, fördelat på olika processor-kärnor.
    chunk_storlek = iterationer_dummy // antal_cores
    chunk_ranges = [(i * chunk_storlek, (i + 1) * chunk_storlek) for i in range(antal_cores)]
    chunk_ranges[-1] = (chunk_ranges[-1][0], iterationer_dummy)

    # Inbakade argument för funktionen.
    args_packed = [(start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
                    slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi,
                    radionuklid_intensitet, radionuklid_sannolikhet) for start, end in chunk_ranges]

    with mp.Pool(antal_cores) as pool:
        partial_results = pool.map(run_MC_multiprocess, args_packed)

    end_time(start)

    #   ----------------------------------------------------------------------
    #   Riktig körning.
    #   ----------------------------------------------------------------------
    print(
        '\n----------------------------------------------------------------------\nACTUAL RUN\n----------------------------------------------------------------------\n')
    start = time.time()

    # Samma multiprocess-kod som ovan.
    chunk_storlek = iterationer_tot // antal_cores
    chunk_ranges = [(i * chunk_storlek, (i + 1) * chunk_storlek) for i in range(antal_cores)]
    chunk_ranges[-1] = (chunk_ranges[-1][0], iterationer_tot)

    # Samma argument, skillnaden är antalet iterationer (fotoner).
    args_packed = [(start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
                    slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi,
                    radionuklid_intensitet, radionuklid_sannolikhet) for start, end in chunk_ranges]

    with mp.Pool(antal_cores) as pool:
        partial_results = pool.map(run_MC_multiprocess, args_packed)

    # Summera resultatmatriserna från respektive process.
    benmärg_matris_deponerad_energi = np.sum(partial_results, axis=0)

    end_time(start)

    # Spara resultatmatrisen i en numpy fil, som sedan går att visualisera i en separat fil.
    np.save('resultat_multiprocess.npy', benmärg_matris_deponerad_energi)
