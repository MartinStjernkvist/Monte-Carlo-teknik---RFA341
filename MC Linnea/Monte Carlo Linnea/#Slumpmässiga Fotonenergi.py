# Slumpmässiga Fotonenergi
import random as rand

rand.random()  # ger slumptal mellan 0 till 1
import numpy as np

np.random.rand()  # ger slumptal mellan 0 till 1)
import matplotlib.pyplot as plt

# Ursprungsenergier från nukliden
Lu177_energ = [208.3, 112.9, 321.3, 249.7, 71.6, 136.7]  # keV

# Intensitet i % för energierna
Lu177_intens = [10.38, 6.2, 0.216, 0.2012, 0.1726, 0.047]

# Gör om intensiteten till tal istället för %
for i in range(len(Lu177_intens)):
    Lu177_intens[i] = Lu177_intens[i] * 10 ** -2

# Skapar en lista med sannolikheten att en viss energi på fotonen skapas
Lu177_sannolik = np.zeros(len(Lu177_intens))

for i in range(len(Lu177_intens)):
    if i == 0:
        Lu177_sannolik[i] = Lu177_intens[i]
    else:
        Lu177_sannolik[i] = Lu177_intens[i] + Lu177_intens[i - 1]

# Tar ett slumpmässigt värde och genererar antalet
Foton_energi = []
Antal_iterationer = 100
for i in range(Antal_iterationer):
    Slump_tal = np.random.rand()  # Slumptal
    if Slump_tal <= Lu177_sannolik[0]:
        Foton_energi.append(Lu177_energ[0])
    elif Slump_tal <= Lu177_sannolik[1]:
        Foton_energi.append(Lu177_energ[1])
    elif Slump_tal <= Lu177_sannolik[2]:
        Foton_energi.append(Lu177_energ[2])
    elif Slump_tal <= Lu177_sannolik[3]:
        Foton_energi.append(Lu177_energ[3])
    elif Slump_tal <= Lu177_sannolik[4]:
        Foton_energi.append(Lu177_energ[4])
    elif Slump_tal <= Lu177_sannolik[5]:
        Foton_energi.append(Lu177_energ[5])
    else:
        continue

# print(Foton_energi)
