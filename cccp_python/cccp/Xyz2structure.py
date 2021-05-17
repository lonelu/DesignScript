import numpy as np
import string

#property
last='  1.00 20.00         ' #Special pdb field	'Occupancy, Temperature factor, Segment identifier, Element symbol'

#parameters
filepath = "XYZ_py.txt"
xyz = np.loadtxt(filepath)
chains = 4
chL = 28

chain_names = list(string.ascii_uppercase[slice(0, chains)])  #"[A, B, C, ...]"

atom_num = 1
res_num = 1

with open('XYZ_py.pdb', 'w') as the_file:
    for i in range(0, chains):
        for j in range(0, chL):     
            x = xyz[i*chL + j, 0]
            y = xyz[i*chL + j, 1]
            z = xyz[i*chL + j, 2]

            # #write N   
            # line = 'ATOM %6d   N  GLY %s %3d    %8.3f%8.3f%8.3f%s%s\n'%(atom_num, chain_names[i], res_num, x, y, z, last, chain_names[i])
            # atom_num = atom_num + 1           
            # the_file.write(line)   

            #write Ca      
            line = 'ATOM %6d  CA  GLY %s %3d    %8.3f%8.3f%8.3f%s%s\n'%(atom_num, chain_names[i], res_num, x, y, z, last, chain_names[i])
            atom_num = atom_num + 1           
            the_file.write(line)  

            # #write C  
            # line = 'ATOM %6d   C  GLY %s %3d    %8.3f%8.3f%8.3f%s%s\n'%(atom_num, chain_names[i], res_num, x, y, z, last, chain_names[i])
            # atom_num = atom_num + 1           
            # the_file.write(line)  

            # #write O
            # line = 'ATOM %6d   O  GLY %s %3d    %8.3f%8.3f%8.3f%s%s\n'%(atom_num, chain_names[i], res_num, x, y, z, last, chain_names[i])
            # atom_num = atom_num + 1           
            # the_file.write(line)    

            res_num = res_num + 1

        the_file.write("TER")  
        

