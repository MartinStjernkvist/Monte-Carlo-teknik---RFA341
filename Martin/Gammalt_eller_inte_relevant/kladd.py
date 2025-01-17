import numpy as np

# print(4.8 * 10**(-3) / 0.998, 'cm')
# print(4.8 * 10**(-3) / 0.998 * 10**(-2) *10**(-6), 'cm')
#
# print(10 * 10**3)
# print(10 * 10**3 / 10**3, 'E3')
#
# print(0.96 **1000)



from scipy.stats import sem

#define dataset
data = [107.2, 113.2, 112.4]

#calculate standard error of the mean
print(sem(data))

print(np.mean(data))
