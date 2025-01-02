from imports import *

'''WIP'''

class attenueringsdata:

    def __init__(self, x,y,z,energi, attenueringsdata_file, anatomi_file):
        self.df_attenueringsdata = pd.read_excel(attenueringsdata_file, index_col=None)

        self.x = x
        self.y = y
        self.z = z
        self.energi = energi
        self.energi_list = self.df_attenueringsdata['E'].to_list()

        self.foto_list = self.df_attenueringsdata['Photoelectric  (cm^2)'].to_list()
        self.compton_list = self.df_attenueringsdata['Compton (cm^2)'].to_list()
        self.rayleigh_list = self.df_attenueringsdata['Rayleigh (cm^2)'].to_list()