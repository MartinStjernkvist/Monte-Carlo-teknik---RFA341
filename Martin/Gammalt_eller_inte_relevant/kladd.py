import numpy as np

# print(4.8 * 10**(-3) / 0.998, 'cm')
# print(4.8 * 10**(-3) / 0.998 * 10**(-2) *10**(-6), 'cm')
#
# print(10 * 10**3)
# print(10 * 10**3 / 10**3, 'E3')
#
# print(0.96 **1000)


from scipy.stats import sem

# define dataset
# data = [107.2, 113.2, 112.4]

data_1 = [3.95, 4.06, 4.26]
data_2 = [8.61, 8.32, 8.61]

# calculate standard error of the mean
print('standardavvikelse')
print(sem(data_1))
print(sem(data_2))

print('medelvärde')
print(np.mean(data_1))
print(np.mean(data_2))


print(4.09 / 4.07)

print(8.51 / 5.22)