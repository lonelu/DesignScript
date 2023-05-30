from pymol import cmd
import glob
import os

workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/'
for _file in glob.glob(os.path.join(workdir, '*.pdb')):
    cmd.load(_file)
seleobjs = cmd.get_names()

target = seleobjs[0] #This cmd also remove the target from seleobjs if in ipdb
mobiles = seleobjs[1:]
print(target)
print(mobiles)
#seleobjs = cmd.get_object_list('(4afl_*)')
for obj in mobiles: cmd.align(obj, target)
for obj in seleobjs: cmd.save(workdir + obj + '_align.pdb', obj)


'''
The cmds below tried to align two pdbs.
'''
from pymol import cmd
import glob
import os

workdir = '/mnt/e/DesignData/_temp/yang/'

for _file in glob.glob(os.path.join(workdir, '*.pdb')):
    cmd.load(_file)
seleobjs = cmd.get_names()


pdb_rmsds = []
for i in range(0, len(seleobjs) - 1):
    rmsd = cmd.align(seleobjs[i], seleobjs[len(seleobjs) - 1])    
    pdb_rmsds.append((seleobjs[i],seleobjs[len(seleobjs) - 1], rmsd))


with open(workdir + '_pymol.txt', 'w') as f:
    for prd in pdb_rmsds:
        f.write(prd[0] + '\t' + prd[1] + '\t' + '\t'.join([str(round(x, 2)) for x in prd[1]]) + '\n')



