'''
To generate the Combs database, one can run the following scripts in order.
a. run_comb.py
b. run_nr.py
c. run_cluster.py
d. run_reps.py

To run the run_comb.py. inpath_rotamer_dir = '/wynton/scratch/nick.polizzi/db_2p8A_0p35rfree_reduced_reps_biounits_pdb_rotamers/'
Those are the result of applying rotalyze to the PDB files of the biological assemblies for the redundancy corrected set of PDB files in the Q-Bits database.  
I have a script for doing that and will send you an example of how it can be used.
You can find the function for doing that (which I adapted from some code Nick sent me) 
at https://github.com/rckormos/combs_ligand_database/blob/main/rotalyze.py, 
and an example of its use at line 467 of https://github.com/rckormos/combs_ligand_database/blob/main/ligand_database.py (edited) 

'''


import os
import pickle
import combs2

#>>> local

probe_dict = {}
_path_to_probe_path = '/mnt/e/DesignData/Combs/vdM_parent_database_20210617/probe_database/pickle_w_val/'
ind = 0
for file in os.listdir(_path_to_probe_path):
    if not '.pkl' in file:
        continue
    probe_dict[ind] = _path_to_probe_path + file
    ind += 1
with open('/mnt/e/DesignData/Combs/vdM_parent_database_20210617/probe_database/probe_dict.pkl', 'wb') as f:
    pickle.dump(probe_dict, f)

cg_dict = combs2.parse.comb.cg_dicts['coo']
_inpath_prody_dir = '/mnt/e/DesignData/Combs/vdM_parent_database_20210617/db_2p8A_0p35rfree_reduced_reps_biounits_prody/'
_inpath_rotamer_dir = '/mnt/e/DesignData/Combs/vdM_parent_database_20210617/rota/'
_outdir = '/mnt/e/DesignData/Combs/vdM_parent_database_20210617/_test_comb/'
_path_to_probe_paths = '/mnt/e/DesignData/Combs/vdM_parent_database_20210617/probe_database/probe_dict.pkl'
ind_first = 0
his_option = None

combs2.parse.comb.run_comb(cg_dict, _inpath_prody_dir, _inpath_rotamer_dir, _outdir, _path_to_probe_paths, ind_first, his_option)



#>>> wynton run_comb
import os
import pickle
import combs2
import multiprocessing as mp
workdir = '/wynton/home/degradolab/lonelu/DesignData/Database/vdM_parent/'

#>>> One time run code to prepare probe_dict.pkl file
# probe_dict = {}
# _path_to_probe_path = workdir + '/pickle_w_val/'
# ind = 0
# for file in os.listdir(_path_to_probe_path):
#     if not '.pkl' in file:
#         continue
#     probe_dict[ind] = _path_to_probe_path + file
#     ind += 1
# with open(workdir + 'probe_dict.pkl', 'wb') as f:
#     pickle.dump(probe_dict, f)

cg_dict = combs2.parse.comb.cg_dicts['coo']
_inpath_prody_dir = workdir + '/db_2p8A_0p35rfree_reduced_reps_biounits_prody/'
_inpath_rotamer_dir = workdir + 'db_2p8A_0p35rfree_reduced_reps_biounits_rot/'
_outdir = workdir + '/coo/'
_path_to_probe_paths = workdir + 'probe_dict.pkl'
ind_first = 0
his_option = None

pool = mp.Pool(48)
[pool.apply_async(combs2.parse.comb.run_comb, args=(cg_dict, _inpath_prody_dir, _inpath_rotamer_dir, _outdir, _path_to_probe_paths, ind_first, his_option)) for ind_first in range(0, 20568)]
pool.close()
pool.join()
#combs2.parse.comb.run_comb(cg_dict, _inpath_prody_dir, _inpath_rotamer_dir, _outdir, _path_to_probe_paths, ind_first, his_option)

# import pandas as pd
# _test = pd.read_parquet(_outdir + '5Y0D_biomol_1_A_D.parquet.gzip')
# _test.to_csv(_outdir + '_test.csv')


#>>> wynton run_nr

import os
import pickle
import combs2

workdir = '/wynton/home/degradolab/lonelu/DesignData/Database/vdM_parent/'

# probe_dict = []
# _path_to_probe_path = workdir + 'pickle_w_val/'
# for file in sorted(os.listdir(_path_to_probe_path))[0:10]:
#     if not '.pkl' in file:
#         continue
#     probe_dict.append(_path_to_probe_path + file) 

# with open(workdir + 'probe_path_list2.pkl', 'wb') as f:
#     pickle.dump(probe_dict, f)


cg = 'coo'
cg_dict = combs2.parse.comb.cg_dicts['coo']

_path_to_probe_paths = workdir + 'probe_path_list2.pkl'
_outdir = workdir + '/_test_nr/'
_path_to_comb_output = workdir + cg + '/'

combs2.parse.nr.run_nr(cg_dict, _path_to_probe_paths, _outdir, _path_to_comb_output)

import pandas as pd
_test = pd.read_parquet(_outdir + 'ALA.parquet.gzip')
_test.to_csv(_outdir + '_test.csv')

#>>> wynton run_cluster
cg = 'coo'
cg_dict = combs2.parse.comb.cg_dicts['coo']
_outdir = workdir + '/_test_comb_cluster/'
_path_to_comb_nr_output = _outdir + cg + '_nr/'
ind = 0
rc = 0.65
cghnp = cgh_name_path = None
min_clu_size = 2
combs2.parse.cluster.run_cluster(cg_dict, _path_to_comb_nr_output, _outdir, ind,
                                    rmsd_cutoff=rc, path_to_cg_H_names=cghnp, min_cluster_size=min_clu_size)


