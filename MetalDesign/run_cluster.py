import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/DesignScript/MetalDesign')
from ligand_database import *

'''
# Download metal containing pdbs.

workdir = "/mnt/e/DesignData/ligands/CA_rcsb/"

filename='all_rcsb.txt'
organize_rcsb_file(workdir)
download_pdb(workdir, filename, resolution = 2.5)
# cd /mnt/e/DesignData/ligands/FE_rcsb/
# gunzip *.gz

'''

workdir = "/mnt/e/DesignData/ligands/all_metal/"

# According to the prody atom flag (http://prody.csb.pitt.edu/manual/reference/atomic/flags.html#flags), NI MN ZN CO CU MG FE CA are not all flag as ion. 
# Note that calcium is read as CA, which is same as alpha carbon in prody selection. 
# So to select calcium, we need to add ion before it.

metal_sel = 'ion or name NI MN ZN CO CU MG FE' 

align_sel_backbone = 'name C CA N O NI MN ZN CO CU MG FE or ion'

aa = 'resname HIS'

'''
# Extract metal core.
# All rcsb pids are extracted by each NI core with sequence (9 aa)

#pdbs = get_all_pbd_prody(workdir + 'all_rcsb/')
#cores = extract_all_core_seq(pdbs, metal_sel, extend = 4)

cores = extract_all_core_seq_from_path(workdir + 'all_rcsb/', metal_sel, extend = 4)

superimpose_core_and_writepdb(cores, cores[0], metal_sel, workdir + '1_seq_cores/')
'''

'''
# cluster based on core sequence. write summary and extract represent (the first on in each cluster)

pdbs = get_all_pbd_prody(workdir + '1_seq_cores/')

clusters = reduce_dup(pdbs, metal_sel)

write_dup_summary(workdir, pdbs, clusters)

extract_rep_and_writepdb(pdbs, clusters, metal_sel, workdir + '2_seq_cores_reps/')

'''

'''
# Futher extract core (metal with contact aa) with aas from cores_seq for manual check.

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

cores = extract_all_core(pdbs, metal_sel)

writepdb(cores, workdir + '3_cores_reps/')

'''

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

# Align his core

aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS')
writepdb(aa_cores, workdir + '4_his_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '4_his_cores_reps/')
run_cluster(_pdbs, workdir, '5_his_core_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 5, align_sel = align_sel_backbone)

# Align his core + phipsi.

aa_ps_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', consider_phipsi = True)
writepdb(aa_ps_cores, workdir + '6_his_ps_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '6_his_ps_cores_reps/')
run_cluster(_pdbs, workdir, '7_his_core_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 9, align_sel = align_sel_backbone)

# Align his core + sidechain

_pdbs = get_all_pbd_prody(workdir + '4_his_cores_reps/')
run_cluster(_pdbs, workdir, '8_his_sc_core_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 11, align_sel = 'heavy')

# Align his core + sidechain + phipsi

_pdbs = get_all_pbd_prody(workdir + '6_his_ps_cores_reps/')
run_cluster(_pdbs, workdir, '9_his_sc_ps_core_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 15, align_sel = 'heavy')

# Align his core + 2 AA 

aa_2aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', extention= 2)
writepdb(aa_2aa_cores, workdir + '10_his_2aa_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '10_his_2aa_cores_reps/')
run_cluster(_pdbs, workdir, '11_his_2aa_cores_cluster/', rmsd = 0.5, metal_sel = metal_sel, len_sel = 21, align_sel = align_sel_backbone)

# Align 2 his core 

aa_2aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', extention= 3, extract2aa=True)
writepdb(aa_2aa_cores, workdir + '12_2his_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '12_2his_cores_reps/')
for i in range(5, 10):
    subworkdir = '13_2his_cores_cluster_' + str(i) + '/'
    run_cluster(_pdbs, workdir, subworkdir, rmsd = 0.5, metal_sel = metal_sel, len_sel = i*4 + 1, align_sel = align_sel_backbone)

# Align atom core

#pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

atom_cores = extract_all_atom_core(pdbs, metal_sel, tag = '_atom')
writepdb(atom_cores, workdir + '14_atom_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '14_atom_cores_reps/')
for i in range(3, 6):
    subworkdir = '14_atom_cores_clu_' + str(i) + '/'
    run_cluster(_pdbs, workdir, subworkdir, rmsd = 0.2, metal_sel = metal_sel, len_sel = i+1, align_sel = 'all')