import os
import prody as pr

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2_coil/_z_fix_c2/'

for x in os.listdir(workdir):
    if not '.pdb' in x:
        continue 
    pdb = pr.parsePDB(workdir + x)
    with open(workdir + x[0:-4] + '.txt', 'w') as f:
        print(x)
        for x in pdb.select('name CA'):
            f.write('\t'.join([str(c) for c in x.getCoords()]) + '\n')
    