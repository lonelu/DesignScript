import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/DesignScript/MetalDesign')
from ligand_database import *

metal_sel = 'name NI or name MN'

workdir = "/mnt/e/DesignData/ligands/NI_MN/"


# Futher extract core with aas from cores_seq for manual check.

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

cores = extract_all_core(pdbs, metal_sel)

writepdb(cores, workdir + '3_cores_reps/')

# Extract his core # Align his core

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS')

writepdb(aa_cores, workdir + '4_his_cores_reps/')

pdbs = get_all_pbd_prody(workdir + '4_his_cores_reps/')

clu = superimpose_aa_core(pdbs, len_sel = 5, align_sel = 'name C CA N O NI MN')

metal_diffs = print_cluster_pdbs(clu, workdir + '5_his_core_cluster/')

write_metal_diff(workdir + '5_his_core_cluster/', metal_diffs)

clu_infos = get_clu_info_write(workdir + '5_his_core_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '5_his_core_cluster/_score.png')

# Align his core + phipsi.

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

aa_ps_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', consider_phipsi = True)

writepdb(aa_ps_cores, workdir + '6_his_ps_cores_reps/')

pdbs = get_all_pbd_prody(workdir + '6_his_ps_cores_reps/')

clu = superimpose_aa_core(pdbs, len_sel = 9, align_sel = 'name C CA N O NI MN')

metal_diffs = print_cluster_pdbs(clu, workdir + '7_his_core_cluster/')

write_metal_diff(workdir + '7_his_core_cluster/', metal_diffs)

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

metal_diffs = print_cluster_pdbs(clu, workdir + '8_his_sc_core_cluster/')

write_metal_diff(workdir + '8_his_sc_core_cluster/', metal_diffs)

clu_infos = get_clu_info_write(workdir + '8_his_sc_core_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '8_his_sc_core_cluster/_score.png')

# Align his core + sidechain + phipsi

pdbs = get_all_pbd_prody(workdir + '6_his_ps_cores_reps/')

clu = superimpose_aa_core(pdbs, len_sel = 15, align_sel = 'heavy')

metal_diffs = print_cluster_pdbs(clu, workdir + '9_his_sc_ps_core_cluster/')

write_metal_diff(workdir + '9_his_sc_ps_core_cluster/', metal_diffs)

clu_infos = get_clu_info_write(workdir + '9_his_sc_ps_core_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '9_his_sc_ps_core_cluster/_score.png')

# Align his core + 2 AA 

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

aa_2aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', extention= 2)

writepdb(aa_2aa_cores, workdir + '10_his_2aa_cores_reps/')

pdbs = get_all_pbd_prody(workdir + '10_his_2aa_cores_reps/')

clu = superimpose_aa_core(pdbs, len_sel = 21, align_sel = 'name C CA N O NI MN')

metal_diffs = print_cluster_pdbs(clu, workdir + '11_his_2aa_cores_cluster/')

write_metal_diff(workdir + '11_his_2aa_cores_cluster/', metal_diffs)

clu_infos = get_clu_info_write(workdir + '11_his_2aa_cores_cluster/_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'heavy')

plot_clu_info(clu_infos, workdir + '11_his_2aa_cores_cluster/_score.png')

# Align 2 his core 

pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

aa_2aa_cores = extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', extention= 3, extract2aa=True)

writepdb(aa_2aa_cores, workdir + '12_2his_cores_reps/')

_pdbs = get_all_pbd_prody(workdir + '12_2his_cores_reps/')

for i in range(5, 10):
    clu = superimpose_aa_core(_pdbs, len_sel = i*4 + 1, align_sel = 'name C CA N O NI MN')

    if not clu or len(clu.mems) == 0: continue
    subworkdir = workdir + '13_2his_cores_cluster_' + str(i) + '/'
    metal_diffs = print_cluster_pdbs(clu, subworkdir)
    write_metal_diff(subworkdir, metal_diffs)
    clu_infos = get_clu_info_write(subworkdir + '_summary.txt', pdbs, clu, rmsd = 0.5, metal_sel = metal_sel, align_sel = 'name C CA N O NI')

    plot_clu_info(clu_infos, subworkdir + '_score.png')


