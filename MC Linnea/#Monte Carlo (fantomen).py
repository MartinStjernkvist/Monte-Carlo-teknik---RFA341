#Monte Carlo (fantomen)
import os

file='FullPhantom.bin.crdownload'
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
abs_file_path = os.path.join(script_dir, file)

x=256
y=256
z=1200
voxelSize= 0.15 #cm

"""""""""""""""
f=open(abs_file_path,"r")
raw_data=f.read(abs_file_path)
f.close()
print(raw_data)

"""""""""""""""

