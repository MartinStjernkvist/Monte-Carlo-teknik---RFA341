#Övningsuppgifter Monte Carlo
import numpy as np

#1
Slump_radvektor=np.zeros((1,10))
print(Slump_radvektor)
for i in range(10):
    Slump_radvektor[0,i]=np.random.normal(0,0.2)

print(Slump_radvektor)

Radvektor=np.random.normal(10)

print(Radvektor)
#2
