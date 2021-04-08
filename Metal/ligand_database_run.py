import os
import prody as pr
import ligand_database as ld

metal_sel = 'name NI'

workdir = "/mnt/e/DesignData/ligands/NI_rcsb/"

# All rcsb pids are extracted by each NI core with sequence (5 aa)

cores, first = ld.extract_all_core_seq(workdir, metal_sel)

ld.superimpose_core_and_write(cores, first, workdir, metal_sel)

# cluster based on core sequence. write summary and extract represent (the first on in each cluster)

pdbs = get_all_pbd_prody(workdir + 'seq_cores/')

clusters = reduce_dup(pdbs, metal_sel)

write_dup_summary(workdir, pdbs, clusters)

extract_rep_and_writepdb(workdir, pdbs, clusters, 'seq_cores_reps/')

# Futher extract core with aas from cores_seq 

cores = extract_all_core(workdir + 'seq_cores_reps/', metal_sel)

writepdb(cores, workdir + 'cores_reps/')

# Extract his core

pdbs = get_all_pbd_prody(workdir + 'cores_reps/')

aa_cores = extract_all_core_aa(pdbs, aa = 'resname HIS')

writepdb(aa_cores, workdir + 'his_cores_reps/')

# Align his core

pdbs = get_all_pbd_prody(workdir + 'his_cores_reps/')

#[p.getTitle() for p in pdbs if len(p)!=11]

superimpose_aa_core(pdbs, workdir + 'his_core_cluster/')