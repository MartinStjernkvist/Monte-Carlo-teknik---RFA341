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
#Intensitet för energierna
Lu177_intens=[10.38, 6.2, 0.216, 0.2012, 0.1726, 0.047]
At211_intens=[41.78, 0.0039, 0.0011, 0.0004, 0.00004]
#Interpolera 
def func(x,a,b,c,d,e,f):

    return a*x**5+b*x**4+c*x**3+d*x**2+e*x+f
params,cv=curve_fit(func,Lu177_energ,Lu177_intens)


def Funktion(x,a,b,c):
   
    return a*np.exp(-b*x)+c 

#Hitta ett sätt att slumpa ursprungenergin