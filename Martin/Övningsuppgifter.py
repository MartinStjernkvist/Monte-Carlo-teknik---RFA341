import numpy as np
from numpy import random

# 1
random_radvektor = np.random.rand(10)


# print(random_radvektor)


# 2
def exponential_fkn(x):
    return 5 ** (10 * x)


list_input = np.random.rand(1000)
list_input = list_input.tolist()
new_list = list(map(exponential_fkn, list_input))


def normalize(x, x_max):
    return x / x_max


x_max = max(new_list)
new_list_normalised = list(map(lambda x: normalize(x, x_max), new_list))
print(new_list_normalised)

# 3
