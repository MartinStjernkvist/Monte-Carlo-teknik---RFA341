
import json
# start = time.time()
#
# for i in range(10**6):
#     print('hej')

# end_time(start)



def inputs_riktig_körning():
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
        iterationer_tot = eval(input_riktig_skalär_iterationer) * 10 ** (eval(input_riktig_magnitud_iterationer))

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

with (open('inputs_upg1_multiprocess.json', 'r') as f):
    json_object = json.load(f)
    antal_cores = json_object['antal_cores']
    iterationer_dummy = json_object['iterationer_dummy']
    iterationer_tot = json_object['iterationer_tot']
    print('output: ', antal_cores, iterationer_dummy, iterationer_tot)
    f.close()