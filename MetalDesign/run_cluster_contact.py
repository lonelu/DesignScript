
import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/DesignScript/MetalDesign')
from ligand_database import *


workdir = "/mnt/e/DesignData/ligands/ZN_rcsb/"

# According to the prody atom flag (http://prody.csb.pitt.edu/manual/reference/atomic/flags.html#flags), NI MN ZN CO CU MG FE CA are not all flag as ion. 
# Note that calcium is read as CA, which is same as alpha carbon in prody selection. 
# So to select calcium, we need to add ion before it.

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE or ion'


pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')


# Align atom core

atom_cores = extract_all_atom_core(pdbs, metal_sel, tag = '_atom')
writepdb(atom_cores, workdir + '8_atom_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '8_atom_cores_reps/')
for i in range(3, 8):
    subworkdir = '8_atom_cores_clu_' + str(i) + '/'
    run_cluster(_pdbs, workdir, subworkdir, rmsd = 0.2, metal_sel = metal_sel, len_sel = i+1, align_sel = 'all')