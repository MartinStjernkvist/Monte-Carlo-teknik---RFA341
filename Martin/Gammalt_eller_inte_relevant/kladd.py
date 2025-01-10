import numpy as np

vec = np.array([1,1,1])

print(vec[1])

print(type(vec[1]))


energi_MeV_list = [1,2,3,4]
energi_list = [(lambda x: x * 10**6)(x) for x in energi_MeV_list]
print(energi_list)