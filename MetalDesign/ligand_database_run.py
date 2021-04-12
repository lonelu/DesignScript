import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/DesignScript/MetalDesign')
from ligand_database import *

metal_sel = 'name NI'

workdir = "/mnt/e/DesignData/ligands/NI_rcsb/"

# All rcsb pids are extracted by each NI core with sequence (7 aa)

pdbs = get_all_pbd_prody(workdir + 'all_rcsb/')

cores = extract_all_core_seq(pdbs, metal_sel, extend = 3)

superimpose_core_and_writepdb(cores, cores[0], metal_sel, workdir + '1_seq_cores/')

# cluster based on core sequence. write summary and extract represent (the first on in each cluster)

pdbs = get_all_pbd_prody(workdir + '1_seq_cores/')

clusters = reduce_dup(pdbs, metal_sel)

write_dup_summary(workdir, pdbs, clusters)

extract_rep_and_writepdb(pdbs, clusters, metal_sel, workdir + '2_seq_cores_reps/')

# Futher extract core with aas from cores_seq for manual check.

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

cores = extract_all_core(pdbs, metal_sel)

writepdb(cores, workdir + '3_cores_reps/')

# Extract his core # Align his core

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS')

writepdb(aa_cores, workdir + '4_his_cores_reps/')

pdbs = get_all_pbd_prody(workdir + '4_his_cores_reps/')

clu = superimpose_aa_core(pdbs, len_sel = 5, align_sel = 'name C CA N O NI')

print_cluster_pdbs(clu, workdir + '5_his_core_cluster/')

clu_infos = get_clu_info_write(workdir + '5_his_core_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '5_his_core_cluster/_score.png')

# Align his core + phipsi.

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

aa_ps_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', consider_phipsi = True)

writepdb(aa_ps_cores, workdir + '6_his_ps_cores_reps/')

pdbs = get_all_pbd_prody(workdir + '6_his_ps_cores_reps/')

clu = superimpose_aa_core(pdbs, len_sel = 9, align_sel = 'name C CA N O NI')

print_cluster_pdbs(clu, workdir + '7_his_core_cluster/')

clu_infos = get_clu_info_write(workdir + '7_his_core_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '7_his_core_cluster/_score.png')

# Align his core + sidechain

pdbs = get_all_pbd_prody(workdir + '4_his_cores_reps/')

'''
_test = []
for p in pdbs:
    if len(p.select('heavy')) != 11:
        _test.append(p.getTitle())
'''

clu = superimpose_aa_core(pdbs,  len_sel = 11, align_sel = 'heavy')

print_cluster_pdbs(clu, workdir + '8_his_sc_core_cluster/')

clu_infos = get_clu_info_write(workdir + '8_his_sc_core_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '8_his_sc_core_cluster/_score.png')

# Align his core + sidechain + phipsi

pdbs = get_all_pbd_prody(workdir + '6_his_ps_cores_reps/')

clu = superimpose_aa_core(pdbs, len_sel = 15, align_sel = 'heavy')

print_cluster_pdbs(clu, workdir + '9_his_sc_ps_core_cluster/')

clu_infos = get_clu_info_write(workdir + '9_his_sc_ps_core_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '9_his_sc_ps_core_cluster/_score.png')

# Align his core + 2 AA 

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

aa_2aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', extention= 2)

writepdb(aa_2aa_cores, workdir + '10_his_2aa_cores_reps/')

pdbs = get_all_pbd_prody(workdir + '10_his_2aa_cores_reps/')

clu = superimpose_aa_core(pdbs, len_sel = 21, align_sel = 'name C CA N O NI')

print_cluster_pdbs(clu, workdir + '11_his_2aa_cores_cluster/')

clu_infos = get_clu_info_write(workdir + '11_his_2aa_cores_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '11_his_2aa_cores_cluster/_score.png')