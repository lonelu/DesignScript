# Mutation to all ala for topology generation.
from metalprot.basic import prody_ext
import prody as pr

workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC_short/loop5/Rosetta/'
target = pr.parsePDB(workdir + 'bb_prep.pdb')
title = 'mut'
mut = prody_ext.target_to_all_gly_ala(target, title, win_no_mutation = [], aa = 'ALA', keep_no_protein = False)
pr.writePDB(workdir + 'bb_prep_ALA_topo.pdb', mut)


# Run topology 
python /mnt/e/GitHub_Design/DesignScript/topology/topology_usage.py