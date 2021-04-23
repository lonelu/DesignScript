import prody as pr

target_path = '/mnt/e/DesignData/ligands/Design_Sam/3into4_helix_assembly_renum.pdb'

target = pr.parsePDB(target_path)

target.select('protein and name CA').getResindices()
target.select('protein and name CA').getResnums().astype(np.int16)
target.select('protein and name CA').getChids().astype('object')

query_path = '/mnt/e/DesignData/ligands/ZN_rcsb/7_HIS_HIS_cores_cluster_5/0/cluster_0_mem_141_centroid_5keb_ZN_1_HIS_HIS_1.pdb'

query = pr.parsePDB(query_path)

rmsd_cut = 1

win_extract = supperimpose_target_bb(query, target, rmsd_cut)

outdir = '/mnt/e/DesignData/ligands/Design_Sam/test/'

write_query_pdbs(outdir, win_extract)

query_2nd_path = '/mnt/e/DesignData/ligands/ZN_rcsb/4_CYS_cores_cluster/0/cluster_0_mem_1335_centroid_4ype_ZN_3_CYS_1.pdb'

query_2nd = pr.parsePDB(query_2nd_path)

win_extract_2nd = supperimpose_target_bb(query_2nd, target, rmsd_cut)

distance_cut = 1

extracts = metal_distance_extract(target, win_extract, win_extract_2nd, distance_cut)

outdir = '/mnt/e/DesignData/ligands/Design_Sam/test_cores/'

write_cores(outdir, win_extract, win_extract_2nd, extracts)