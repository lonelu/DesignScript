# The surface of a protein is related to the solubility and cystallization. 
# This script is to measure the number of K R D E of a protein and calculate the charge.
import os
import prody as pr

workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC_redesign/redesign_599_surface2/'

pdbs = []

for file in os.listdir(workdir):
    if not '.pdb' in file:
        continue
    pdbs.append(pr.parsePDB(workdir + file))

for pdb in pdbs:
    Kc = len(pdb.select('name CA and resname LYS'))
    Rc = len(pdb.select('name CA and resname ARG'))
    Ec = len(pdb.select('name CA and resname GLU'))
    Dc = len(pdb.select('name CA and resname ASP'))
    charge = Kc + Rc - Ec - Dc
    print(pdb.getTitle() + '\t' + str(Kc) + '\t'+ str(Rc) + '\t' + str(Ec) + '\t' + str(Dc) + '\t' + str(charge)) 