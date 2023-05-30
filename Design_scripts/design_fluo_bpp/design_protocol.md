### install miniconda.

### install conda environment.

### prepare .params file.
# prepare .params file on local.
cd /mnt/e/DesignData/bpp_fluo_comb/44b/
/mnt/d/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n 44b -p 44b --conformers-in-one-file 44b.sdf
/mnt/d/Rosetta/main/source/scripts/python/public/molfile_to_params.py -n dhb -p dhb --conformers-in-one-file dhb.sdf

# prepare .params file on wynton.
/wynton/home/degradolab/lonelu/rosetta_bin_linux_2021.16.61629_bundle/main/source/scripts/python/public/molfile_to_params.py -n dhb -p dhb --conformers-in-one-file dhb.sdf

### run COMBS2


### run Rosetta

# Mutation to all ala for topology generation.

from metalprot.basic import prody_ext
import prody as pr

workdir = '/mnt/e/DesignData/bpp_fluo_comb/fluo/output1_09_f63440_nick_ala/sel/Rosetta/'
target = pr.parsePDB(workdir + 'bb_only.pdb')
title = 'mut'
mut = prody_ext.target_to_all_gly_ala(target, title, win_no_mutation = [], aa = 'ALA', keep_no_protein = False)
pr.writePDB(workdir + 'bb_prep_ALA_topo.pdb', mut)

# Use Rosetta to mutate.
python /mnt/e/GitHub_Design/DesignScript/pyrosetta/pyrosetta_mutation.py

# Run topology 
python /mnt/e/GitHub_Design/DesignScript/topology/topology_usage.py

 
### run OmegaFold

