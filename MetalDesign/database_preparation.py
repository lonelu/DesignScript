import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/DesignScript/MetalDesign')
from ligand_database import *


# Download metal containing pdbs.

workdir = "/mnt/e/DesignData/ligands/CA_rcsb/"

filename='all_rcsb.txt'
organize_rcsb_file(workdir)
download_pdb(workdir, filename, resolution = 2.5)
# cd /mnt/e/DesignData/ligands/FE_rcsb/
# gunzip *.gz



#------------------------------------------------------------------


# Extract metal core.
# All rcsb pids are extracted by each NI core with sequence (9 aa)

workdir = "/mnt/e/DesignData/ligands/ZN_rcsb/"
metal_sel = 'name ZN'

#pdbs = get_all_pbd_prody(workdir + 'all_rcsb/')
#cores = extract_all_core_seq(pdbs, metal_sel, extend = 4)

cores = extract_all_core_seq_from_path(workdir + 'all_rcsb/', metal_sel, extend = 4)

superimpose_core_and_writepdb(cores, cores[0], metal_sel, workdir + '1_seq_cores/')


#------------------------------------------------------------------


# cluster based on core sequence. write summary and extract represent (the first on in each cluster)

pdbs = get_all_pbd_prody(workdir + '1_seq_cores/')

clusters = reduce_dup(pdbs, metal_sel)

write_dup_summary(workdir, pdbs, clusters)

extract_rep_and_writepdb(pdbs, clusters, metal_sel, workdir + '2_seq_cores_reps/')

