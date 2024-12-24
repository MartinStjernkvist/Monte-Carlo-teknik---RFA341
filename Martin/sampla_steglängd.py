from imports import *

def invers_funktion(x, mu):
    # x borde i det här fallet vara rand(0,1)
    return -np.log(x) / mu

def medelvägslängd(mu):
    medelvägslängd = invers_funktion(np.random.rand(), mu)

# FELAKTIGT:
# def medelvägslängd_lista(iterationer, mu):
#     medelvägslängd_lista = np.zeros(iterationer)
#     for i in range(iterationer):
#         medelvägslängd_lista = invers_funktion(np.random.rand(), mu)



#   ----------------------------------------------------------------------
#   TEST, KOMMENTERA BORT
#   ----------------------------------------------------------------------

"""
mu_max = 1 # placeholder, byt till maximala mu som finns för fotonenergierna
iterationer = 100
R = np.zeros(iterationer)
steglängd_list = np.zeros(iterationer)

for i in range(iterationer):
    R[i] = np.random.rand()
    steglängd_list[i] = invers_funktion(R[i], mu_max)

print(steglängd_list)

x_list = np.linspace(0,1,1000)


plt.plot(x_list, invers_funktion(x_list, mu_max))
plt.show()
"""