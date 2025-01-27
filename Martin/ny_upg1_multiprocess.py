from imports import *

from upg1_sampla_energi_start import energi_start
from upg1_matriser import slicad_fantom_matris, slicad_njure_matris, slicad_benmärg_matris
from upg1_sampla_position_start import position_start
from upg1_sampla_riktning_och_steg_start import riktning_uniform, steg
from upg1_sampla_steglängd import medelvägslängd
from upg1_sampla_växelverkan import växelverkan
from upg1_attenueringsdata import attenueringsdata
from upg1_sampla_compton import compton_vinkel_och_energiförlust
from upg1_sampla_foto_vxv import foto_vxv
from upg1_bestäm_om_attenuerad import bestäm_om_attenuerad
from upg12_förflyttning import förflyttning
from upg1_bestäm_om_vxv import bestäm_om_vxv_sker

# from upg1_steg_transformation import ny_steg_transformera_koordinatsystem_3d
from upg12_steg_transformation import transformera_koordinatsystem
from upg12_rotation_matris import rotations_matris


def steg_tills_vxv(x, y, z, theta, phi, steglängd, R, foton_energi, mu_max, x_size, y_size, z_size, df_attenueringsdata,
                   df_anatomidefinitioner, utanför_fantom):
    x_start, y_start, z_start = x, y, z
    vxv_sker = False

    while vxv_sker == False:
        # Sampla medelvägslängden från inverstransformerad attenueringsfunktion.
        steglängd_ny = medelvägslängd(mu_max)
        # print(f'steglängd: {steglängd}')

        # Gå steget till ny position från startpositionen i startriktningen.
        # Eftersom voxlarna har diskreta positioner måste avrundning till närmaste heltal göras.

        dx, dy, dz, R = transformera_koordinatsystem(steglängd, phi,
                                                     theta,
                                                     steglängd_ny,
                                                     0,
                                                     0, R)

        x_ny, y_ny, z_ny, x_round, y_round, z_round = förflyttning(x_start, y_start, z_start, dx, dy, dz,
                                                                   voxel_sidlängd)

        x_start, y_start, z_start = x_ny, y_ny, z_ny
        steglängd = steglängd_ny

        #   -----------------------------------
        #   Om foton hamnar utanför fantommatrisen
        #   -> kasta ut foton ur loopen.
        #   Ekvationen förekommer nedan efter varje nytt steg tas.
        #   -----------------------------------
        attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size,
                                                          utanför_fantom, slicad_fantom_matris, foton_energi)

        if attenuerad == 0:
            voxelvärde = slicad_fantom_matris[x_round, y_round, z_round]

            vxv_sker = bestäm_om_vxv_sker(voxelvärde, foton_energi, mu_max, df_attenueringsdata,
                                          df_anatomidefinitioner)
        else:
            # utanför matris
            vxv_sker = True

    """
    while vxv_sker == False:
        # Sampla medelvägslängden från inverstransformerad attenueringsfunktion.
        steglängd = medelvägslängd(mu_max)
        # print(f'steglängd: {steglängd}')

        # Gå steget till ny position från startpositionen i startriktningen.
        # Eftersom voxlarna har diskreta positioner måste avrundning till närmaste heltal göras.
        dx, dy, dz = steg(theta, phi, steglängd)
        x_ny, y_ny, z_ny, x_round, y_round, z_round = förflyttning(x_start, y_start, z_start, dx, dy, dz, voxel_sidlängd)

        x_start, y_start, z_start = x_ny, y_ny, z_ny

        #   -----------------------------------
        #   Om foton hamnar utanför fantommatrisen
        #   -> kasta ut foton ur loopen.
        #   Ekvationen förekommer nedan efter varje nytt steg tas.
        #   -----------------------------------
        attenuerad, utanför_fantom = bestäm_om_attenuerad(x_round, y_round, z_round, x_size, y_size, z_size,
                                                          utanför_fantom, slicad_fantom_matris, foton_energi)

        if attenuerad == 0:
            voxelvärde = slicad_fantom_matris[x_round, y_round, z_round]

            vxv_sker = bestäm_om_vxv_sker(voxelvärde, foton_energi, mu_max, df_attenueringsdata,
                                          df_anatomidefinitioner)
        else:
            # utanför matris
            vxv_sker = True
    """

    # print(f'attenuerad: {attenuerad}')
    return attenuerad, x_ny, y_ny, z_ny, x_round, y_round, z_round, steglängd


def run_MC_multiprocess(args):
    (start, end, tvärsnitt_file, attenueringsdata_file, anatomidefinitioner_file, slicad_fantom_matris,
     slicad_njure_matris, slicad_benmärg_matris, voxel_sidlängd, radionuklid_energi, radionuklid_intensitet,
     radionuklid_sannolikhet) = args

    #   -----------------------------------
    #   Läs in data.
    #   -----------------------------------
    df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)
    df_anatomidefinitioner = pd.read_excel(anatomidefinitioner_file, index_col=None)
    df_tvärsnitt = pd.read_excel(tvärsnitt_file, index_col=None)

    # Matrisdimensionerna, viktigt för att bedöma ifall fotonen är innanför matrisen eller inte.
    x_size, y_size, z_size = slicad_fantom_matris.shape

    # Skapa tom matris, för att sedan fylla på med energideponering i voxlarna.
    benmärg_matris_deponerad_energi = np.zeros((x_size, y_size, z_size))

    #   -----------------------------------
    #   Håll reda på vad som händer när koden körs.
    #   -----------------------------------
    utanför_fantom = 0
    vxv_foto = 0
    träff_benmärg = 0

    #   -----------------------------------
    #   Loopar igenom alla sönderfall (iterationer).
    #   -----------------------------------
    for i in range(start, end):

        # Initiera loopen med att attenuerad = 0.
        attenuerad = 0
        # print(i)

        # Start: sampla position, riktning och energi
        foton_energi = energi_start(radionuklid_energi, radionuklid_sannolikhet)
        x_start, y_start, z_start = position_start(slicad_njure_matris)
        theta, phi = riktning_uniform()
        R = rotations_matris(phi, theta)

        voxelvärde = slicad_fantom_matris[x_start, y_start, z_start]
        instans = attenueringsdata(voxelvärde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
        mu_max = instans.mu_max()

        x, y, z = x_start, y_start, z_start
        steglängd = 0

        #   -----------------------------------
        #   Loopa under tiden som fotonen inte attenuerats
        #   och
        #   fortfarande är i matrisen.
        #   -----------------------------------
        while attenuerad == 0:

            attenuerad, x, y, z, x_round, y_round, z_round, steglängd = steg_tills_vxv(x, y, z,
                                                                                       theta, phi, steglängd, R,
                                                                                       foton_energi,
                                                                                       mu_max,
                                                                                       x_size,
                                                                                       y_size,
                                                                                       z_size,
                                                                                       df_attenueringsdata,
                                                                                       df_anatomidefinitioner,
                                                                                       utanför_fantom)

            if attenuerad == 1:
                break

            # Identifiera vilken voxel fotonen befinner sig i.
            voxelvärde = slicad_fantom_matris[x_round, y_round, z_round]
            instans = attenueringsdata(voxelvärde, foton_energi, df_attenueringsdata, df_anatomidefinitioner)
            mu_max = instans.mu_max()

            # Bestäm vilken typ av växelverkan som sker vid nya positionen.
            instans = växelverkan(foton_energi, df_tvärsnitt)
            vxv = instans.bestäm_växelverkan()
            # print(f'energi: {foton_energi * 10 ** (-3)} keV, vxv: {vxv}')

            #   -----------------------------------
            #   Fotoabsorption.
            #   -----------------------------------
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
                        f'foto: {energi_deponering:.0f} eV i voxel [{x_round, y_round, z_round}]')

                # Om fluorescens sker -> följ ny foton.
                if attenuerad == 0:
                    # Sampla spridningsvinklar (uniformt samplade).
                    theta_foto, phi_foto = riktning_uniform()

                    # Ingångsvärden till koordinat-transformeringen som behöver genomföras ifall växelverkan = compton eller rayleigh.
                    theta, phi = theta_foto, phi_foto
                    R = rotations_matris(phi, theta)

            #   -----------------------------------
            #   Comptonspridning.
            #   -----------------------------------
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
                        f'compton: {energideponering_compton:.0f} eV i voxel [{x_round, y_round, z_round}]')

                steglängd_compton = 0.001
                # Koordinattransformation, eftersom spridningsvinklarna för Comptonspridning inte kan samplas uniformt.
                dx_compton, dy_compton, dz_compton, R = transformera_koordinatsystem(steglängd, phi,
                                                                                     theta,
                                                                                     steglängd_compton,
                                                                                     phi_compton,
                                                                                     theta_compton, R)

                x, y, z, x_round, y_round, z_round = förflyttning(x, y, z, dx_compton, dy_compton, dz_compton,
                                                                  voxel_sidlängd)

                # Ingångsvärden till koordinat-transformeringen (om nästa växelverkan är Comptonspridning eller Rayleighspridning).
                theta, phi = theta_compton, phi_compton
                R = rotations_matris(phi, theta)

            #   -----------------------------------
            #   Rayleighspridning.
            #   -----------------------------------
            elif vxv == 'rayleigh':

                # Sampla spridningsvinklar och steglängd.
                theta_rayleigh = np.arccos(-1 + 2 * np.random.rand())  # ERSÄTT MED THOMSON TVÄRSNITT
                phi_rayleigh = 2 * pi * np.random.rand()

                steglängd_rayleigh = 0.001
                # Koordinattransformation, eftersom spridningsvinklarna för Comptonspridning inte kan samplas uniformt.
                dx_rayleigh, dy_rayleigh, dz_rayleigh, R = transformera_koordinatsystem(steglängd, phi,
                                                                                        theta,
                                                                                        steglängd_rayleigh,
                                                                                        phi_rayleigh,
                                                                                        theta_rayleigh, R)

                x, y, z, x_round, y_round, z_round = förflyttning(x, y, z, dx_rayleigh, dy_rayleigh, dz_rayleigh,
                                                                  voxel_sidlängd)

                # Ingångsvärden till koordinat-transformeringen (om nästa växelverkan är Comptonspridning eller Rayleighspridning).
                theta, phi = theta_rayleigh, phi_rayleigh
                R = rotations_matris(phi, theta)

    # Några print statements, för att hålla reda på vad som händer när koden kör.
    print(f'max värdet av matrisen: {np.max(benmärg_matris_deponerad_energi)}')
    print(f'utanför: {utanför_fantom}')
    print(f'foto: {vxv_foto}')
    print(f'träffar: {träff_benmärg}')

    return benmärg_matris_deponerad_energi


def inputs_riktig_körning():
    print('\n-----------------------------------\nDags att ange parametrar.\n-----------------------------------\n')
    print('Standard eller inte?')
    input_standard = input('\nOm standard: s, Annars: vad som helst: ')

    if input_standard == 's':
        antal_cores = 8
        iterationer_dummy = 10 ** 3
        iterationer_tot = 10 ** 5

    else:
        print(
            '\n-----------------------------------\nVIKTIGT:\n-----------------------------------\nAnge antal processor kärnor')
        input_antal_cores = input('Antal kärnor: ')

        if eval(input_antal_cores) > 8:
            antal_cores = 1
        else:
            antal_cores = eval(input_antal_cores)

        print(
            '\n-----------------------------------\nDUMMY:\n-----------------------------------\nAnge magnitud: ex 3 -> 10^3 iterationer')
        input_dummy_magnitud_iterationer = input('Magnitud: ')

        if eval(input_dummy_magnitud_iterationer) > 4:
            iterationer_dummy = 10 ** 3
        else:
            iterationer_dummy = 10 ** (eval(input_dummy_magnitud_iterationer))

        print(
            '\n-----------------------------------\nRIKTIG:\n-----------------------------------\nAnge skalär och magnitud: ex 5 och 5 -> 5 * 10^5 iterationer')
        input_riktig_skalär_iterationer = input('Skalär: ')
        input_riktig_magnitud_iterationer = input('Magnitud: ')

        if eval(input_riktig_magnitud_iterationer) >= 8:
            iterationer_tot = 10 ** 3
        else:
            iterationer_tot = eval(input_riktig_skalär_iterationer) * 10 ** (
                eval(input_riktig_magnitud_iterationer))

    print('antal_cores, iterationer_dummy, iterationer_tot: ', antal_cores, iterationer_dummy, iterationer_tot)
    return antal_cores, iterationer_dummy, iterationer_tot


def spara_resultat(matris, json_object):
    print(
        '\n-----------------------------------\nDags att sparan resultaten i en matris.\n-----------------------------------\n')
    print(
        'Om du anger namn: skriver du "text" utan citattecken kommer filerna \n[resultat_text.npy] och [inputs_text.json] \nskapas.')
    input_spara_resultat = input('\nVar vill du spara matrisen? Om standard: s, Annars: ange namn: ')

    if input_spara_resultat == 's':
        fil_namn_npy = 'resultat_multiprocess.npy'
        # Spara resultatmatrisen i en numpy fil.
        # Resultatmatrisen går att visualisera i en separat fil.
        np.save(fil_namn_npy, matris)

        fil_namn_json = 'inputs_resultat_multiprocess.json'
        with open(fil_namn_json, 'w') as f:
            f.write(json_object)
            f.close()
    else:
        fil_namn_npy = 'resultat_' + input_spara_resultat + '.npy'
        np.save(fil_namn_npy, matris)

        fil_namn_json = 'inputs_' + input_spara_resultat + '.json'
        with open(fil_namn_json, 'w') as f:
            f.write(json_object)
            f.close()

    return print(f'Resultat sparade i filerna [{fil_namn_npy}] och [{fil_namn_json}]')


if __name__ == "__main__":
    radionuklid_energi = Lu177_energi
    radionuklid_intensitet = Lu177_intensitet
    radionuklid_sannolikhet = Lu177_sannolikhet

    antal_cores, iterationer_dummy, iterationer_tot = inputs_riktig_körning()

    dictionary = {
        "antal_cores": antal_cores,
        "iterationer_dummy": iterationer_dummy,
        "iterationer_tot": iterationer_tot
    }

    json_object = json.dumps(dictionary)

    #   -----------------------------------
    #   Dummy run.
    #   För att snabba på den riktiga körningen av koden.
    #   -----------------------------------
    print(
        '\n-----------------------------------\nDUMMY RUN\n-----------------------------------\n')
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

    #   -----------------------------------
    #   Riktig körning.
    #   -----------------------------------
    print(
        '\n-----------------------------------\nACTUAL RUN\n-----------------------------------\n')
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

    spara_resultat(benmärg_matris_deponerad_energi, json_object)
