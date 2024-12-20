#Monte Carlo Slumpa tal
import random as rand
rand.random()# ger slumptal mellan 0 till 1
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#Ursprungsenergier från nuklider
Lu177_energ=[208.3, 112.9, 321.3, 249.7, 71.6, 136.7] #keV
Y90_energ=2278.7 #keV
At211_energ=[5869, 5211.9,5140.3, 4993.4,4895.4] #keV
#Intensitet i % för energierna
Lu177_intens=[10.38, 6.2, 0.216, 0.2012, 0.1726, 0.047]
At211_intens=[41.78, 0.0039, 0.0011, 0.0004, 0.00004]

#Gör om till tal och tar totala sannolikheten för varje energi från en tom matris
Lu177_sannolik=np.zeros(len(Lu177_intens))
At211_sannolik=np.zeros(len(At211_intens))

for i in range(len(Lu177_intens)):
    Lu177_intens[i]=Lu177_intens[i]*10**-2
    At211_intens[i]=At211_intens[i]*10**-2


for i in range(len(Lu177_intens)):
    if i==0:
        Lu177_sannolik[i]=Lu177_intens[i]
        At211_sannolik[i]=At211_intens[i]
    else:
        Lu177_sannolik[i]= Lu177_intens[i]+Lu177_intens[i-1]
        At211_sannolik[i]= At211_intens[i]+At211_intens[i-1]
#Hitta ett sätt att slumpa ursprungenergin