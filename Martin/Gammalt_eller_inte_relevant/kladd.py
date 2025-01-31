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

data_1 = [4.16, 4.15, 4.16]
data_2 = [8.56, 8.63, 8.68]

# calculate standard error of the mean
print('standardavvikelse')
print('yt: ', sem(data_1))
print('homogen: ', sem(data_2))

print('medelvärde')
print('yt: ',np.mean(data_1))
print('homogen: ', np.mean(data_2))


print(4.09 / 4.07)
print(4.16 / 4.07)
print()
print(8.51 / 5.22)
print(8.62 / 5.22)