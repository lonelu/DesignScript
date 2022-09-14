'''
using Rian's code from https://github.com/rckormos/combs_ligand_database

'''

import os 
import sys
from rotalyze import parse_rotalyze
import multiprocessing as mp

_rotalyze = '/wynton/home/degradolab/lonelu/software/molprobity/build/bin/phenix.rotalyze'

workdir = '/wynton/home/degradolab/lonelu/DesignData/Database/vdM_parent/'

indir = workdir + 'db_2p8A_0p35rfree_reduced_reps_biounits_pdb/'
rotalyze_outdir = workdir + 'db_2p8A_0p35rfree_reduced_reps_biounits_rot/'
os.makedirs(rotalyze_outdir, exist_ok = True)

pdbs = []
for pdb in os.listdir(indir):
    if not pdb.endswith('.pdb'):
        continue
    pdbs.append(pdb)

pool = mp.Pool(48)
[pool.apply_async(parse_rotalyze, args=(indir + pdb, _rotalyze, rotalyze_outdir)) for pdb in pdbs]
pool.close()
pool.join()
#parse_rotalyze(indir + pdb, _rotalyze, rotalyze_outdir + '/' + pdb[1:3])