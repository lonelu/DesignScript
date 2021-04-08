# work with ProtCID database
'''
import os
import ligand_database as ld

workdir = "/mnt/e/DesignData/ligands/NI/"
metal_sel = 'name NI'

all_pdb_infos = ld.extract_pml_pdb_info(workdir, metal_sel)

#ld.write_pdb_info(workdir + 'protcid_info.txt', all_pdb_infos)

clean_pdb_infos = [info for info in all_pdb_infos if info.selection == 1]

copy_pfam_subfolder(workdir, 'pfam2pdbs_sel/', clean_pdb_infos)


'''
# Compare with rcsb database
'''
pdb_names, pdb_not_in_protcid = ld.compare_rcsb_file(workdir, all_pdb_infos)

workdir_rcsb = '/mnt/e/DesignData/ligands/NI_rcsb/'

pdb_info_not_protcid = ld.extract_rcsb_pdb_info(workdir_rcsb, metal_sel, True, pdb_not_in_protcid)

ld.write_pdb_info(workdir_rcsb + 'pdb_info_not_protcid.txt', pdb_info_not_protcid)

'''