# Mutation to all ala for topology generation.
from metalprot.basic import prody_ext
import prody as pr

'''
conda activate env_conda
change the workdir and backbone name.
copy and run in ipython
>ipython
'''

workdir = '/mnt/e/DesignData/bpp_fluo_comb/fluo/output1_09_f63440_nick_ala/sel/Rosetta/'
target = pr.parsePDB(workdir + 'bb_only.pdb')
title = 'mut'
mut = prody_ext.target_to_all_gly_ala(target, title, win_no_mutation = [], aa = 'ALA', keep_no_protein = False)
pr.writePDB(workdir + 'bb_prep_ALA_topo.pdb', mut)
