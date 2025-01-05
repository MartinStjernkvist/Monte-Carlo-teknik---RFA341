from imports import *

# start = time.time()
#
# for i in range(10**6):
#     print('hej')

# end_time(start)



df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)
energi_list = df_attenueringsdata['E'].to_list()

big_list_names = ['water', 'muscle', 'lung', 'dry_spine', 'dry_rib', 'blood', 'heart', 'kidney', 'liver',
                          'lymph', 'pancreas', 'intestine',
                          'cartilage', 'brain', 'spleen', 'air', 'breast_mammary', 'skin', 'eye_lens', 'red_marrow',
                          'yellow_marrow', 'thyroid',
                          'bladder']
energi = 50000

energi_list = np.array(energi_list)
diff = np.abs(energi_list - energi)
closest_indices = np.argsort(diff)[:2]
print(closest_indices)


mu_array = np.zeros(len(big_list_names), dtype=object)
mu_values = np.zeros(len(big_list_names), dtype=float)


for i in range(len(big_list_names)):
    mu_array[i] = np.array(df_attenueringsdata[big_list_names[i]][closest_indices].to_list())
    print(mu_array)

mu_max_close_values = [np.max(arr) for arr in mu_array]
mu_max_close = max(mu_max_close_values)
mu_max_close_index = mu_max_close_values.index(mu_max_close)


print(mu_max_close_index)

print(mu_array[mu_max_close_index])

print(mu_max_close)