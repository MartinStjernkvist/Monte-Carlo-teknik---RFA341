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

# 1     &127,2 \\
# 2& 123,5\\
# 3& 127,6\\

"""
ELEKTRONER
"""
print('BETA')
# data_1 = [4.16, 4.15, 4.16]
# data_2 = [8.56, 8.63, 8.68]
data_1 = [4.178, 4.215, 4.194]
data_2 = [8.648, 8.532, 8.549]

# calculate standard error of the mean
print('standardavvikelse')
print('yt: ', sem(data_1))
print('homogen: ', sem(data_2))

print('medelvärde')
print('yt: ', np.mean(data_1))
print('homogen: ', np.mean(data_2))

print(4.20 / 4.07)
# print(4.16 / 4.07)
# print()
# print(8.51 / 5.22)
print(8.58 / 5.22)

"""
ALFA
"""
print('\nALFA')
data_1 = [3.972, 3.976, 3.982]
data_2 = [22.365, 22.377, 22.269]

# calculate standard error of the mean
print('standardavvikelse')
print('yt: ', sem(data_1))
print('homogen: ', sem(data_2))

print('medelvärde')
print('yt: ', np.mean(data_1))
print('homogen: ', np.mean(data_2))

print(22.3 / 9.18)
print(3.97 / 1.66)
