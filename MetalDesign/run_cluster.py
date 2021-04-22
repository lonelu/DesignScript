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

aa = 'resname CYS' #Change len_sel according to aa in the following code.

aa_name = 'CYS' 

len_sel = 5

len_sel_sc = len_sel + 2 #Change len_sel according to aa sidechain number. His:6, Glu: 5, Asp: 4, Cyc: 2

len_sel_ps = len_sel + 4 

len_sel_ps_sc = len_sel_ps + 2 #Change len_sel according to aa sidechain number. His:6, Glu: 5, Asp: 4, Cyc: 2



'''
# Extract aa core (metal with contact aa) with aas from cores_seq for manual check.

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

cores = extract_all_core(pdbs, metal_sel)

writepdb(cores, workdir + '3_aa_cores_reps/')

'''

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

# Align aa core

aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = aa)
writepdb(aa_cores, workdir + '4_' + aa_name + '_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '4_' + aa_name + '_cores_reps/')
run_cluster(_pdbs, workdir, '4_' + aa_name + '_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel, align_sel = align_sel_backbone)

# Align aa core + sidechain

run_cluster(_pdbs, workdir, '4_' + aa_name + '_sc_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_sc, align_sel = 'heavy') 

# Align aa core + phipsi.

aa_ps_cores = extract_all_core_aa(pdbs, metal_sel, aa = aa, consider_phipsi = True)
writepdb(aa_ps_cores, workdir + '5_' + aa_name + '_ps_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '5_' + aa_name + '_ps_cores_reps/')
run_cluster(_pdbs, workdir, '5_' + aa_name + '_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_ps, align_sel = align_sel_backbone)

# Align aa core + sidechain + phipsi

run_cluster(_pdbs, workdir, '5_' + aa_name + '_sc_ps_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = len_sel_ps_sc, align_sel = 'heavy') 

# Align his core + 2 AA 

aa_2aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = aa, extention= 2)
writepdb(aa_2aa_cores, workdir + '6_' + aa_name + '_2aa_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '6_' + aa_name + '_2aa_cores_reps/')
run_cluster(_pdbs, workdir, '6_' + aa_name + '_2aa_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 21, align_sel = align_sel_backbone)

# Align 2 his core 

aa_2aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = aa, extention= 3, extract2aa=True)
writepdb(aa_2aa_cores, workdir + '7_2' + aa_name + '_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '7_2' + aa_name + '_cores_reps/')
for i in range(5, 10):
    subworkdir = '7_2' + aa_name + '_cores_cluster_' + str(i) + '/'
    run_cluster(_pdbs, workdir, subworkdir, rmsd = 0.5, metal_sel = metal_sel, len_sel = i*4 + 1, align_sel = align_sel_backbone)

# Align atom core

'''
atom_cores = extract_all_atom_core(pdbs, metal_sel, tag = '_atom')
writepdb(atom_cores, workdir + '8_atom_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '8_atom_cores_reps/')
for i in range(3, 6):
    subworkdir = '8_atom_cores_clu_' + str(i) + '/'
    run_cluster(_pdbs, workdir, subworkdir, rmsd = 0.2, metal_sel = metal_sel, len_sel = i+1, align_sel = 'all')
'''