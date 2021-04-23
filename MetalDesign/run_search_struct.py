import os
import sys
import prody as pr
sys.path.append(r'/mnt/e/GitHub_Design/DesignScript/MetalDesign')
from search_struct import *

workdir = '/mnt/e/DesignData/ligands/Design_Sam/'

outdir = workdir + 'output/'

if not os.path.exists(outdir):
    os.mkdir(outdir)

target_path = workdir + '3into4_helix_assembly_renum.pdb'

query_dir = '/mnt/e/DesignData/ligands/ZN_rcsb/'

target = pr.parsePDB(target_path)

#Get query pdbs 

subfolders_with_paths = [f.path for f in os.scandir(query_dir) if f.is_dir() if f.name[0] == '7' and 'cluster' in f.name]

query_pdb_paths = []

for subfolder in subfolders_with_paths:
    paths = extract_centroid_pdb_in_clu(subfolder + '/', file_path = '_summary.txt', score_cut = 0, clu_num_cut = 10)
    query_pdb_paths.extend(paths)

querys = [pr.parsePDB(path) for path in query_pdb_paths]

#Get query_2nd pdbs 

subfolders_with_paths2 = [f.path for f in os.scandir(query_dir) if f.is_dir() if f.name[0] == '6' and 'cluster' in f.name]

query_2nd_pdb_paths = []

for subfolder in subfolders_with_paths2:
    paths = extract_centroid_pdb_in_clu(subfolder + '/', file_path = '/_summary.txt', score_cut = 1, clu_num_cut = 50)
    query_2nd_pdb_paths.extend(paths)

query_2nds = [pr.parsePDB(path) for path in query_2nd_pdb_paths]

# supperimpose
rmsd_cut = 0.5
rmsd_cut_2nd = 0.5
distance_cut = 1
count = 1

for query in querys:
    for query_2nd in query_2nds:
        outpath = outdir + str(count) + '/'
        query_win, query_2nd_win, extracts = build_core(outpath, target, query, query_2nd, rmsd_cut, rmsd_cut_2nd, distance_cut)
        if query_win == 0:
            break
        if extracts and len(extracts) > 0:
            count += 1



    
