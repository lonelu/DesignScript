'''
This is the command I use for the reduce program, where BUILD adds hydrogens, FLIP rotates and flips asn gln his hydrogens, (I think having both BUILD and FLIP are redundant but w/e), and Quiet suppresses verbose output (optional):
reduce -BUILD -FLIP -Quiet oldpdb.pdb > newpdb.PDB

I just realized the gpu computer doesn’t have reduce anymore. Or if it does, idk where the path is. However, it’s in the wynton amber installation: /wynton/home/grabe/shared/amber/amber18/bin/reduce
I don’t know if this is necessary, but you might have to source amber first: source /wynton/home/grabe/shared/amber/amber18/amber.sh
'''

import os 

# Define pre-hydrogen and post-hydrogen database paths
db_dir = '/home/sophia/metals_lei/Database/'
db_H_dir = '/home/sophia/metals_lei/DatabaseH/'

for vdm_dir in os.listdir(db_dir):
    vdm_dirpath = db_dir+vdm_dir+'/'
    # make databaseH dir
    if vdm_dir not in os.listdir(db_H_dir):
        os.mkdir(db_H_dir+vdm_dir)
    for cluster in os.listdir(vdm_dirpath):
        # make sure cluster directory isn't txt or png file
        if cluster.endswith('.txt') or cluster.endswith('png'):
            continue
        # make databaseH dir
        if cluster not in os.listdir(db_H_dir+vdm_dir):
            os.mkdir(db_H_dir+vdm_dir+'/'+cluster)
        
        # run reduce on all pdbs
        clusterpath = vdm_dirpath+cluster
        for pdb in os.listdir(clusterpath):
            pdbpath = clusterpath+'/'+pdb
            os.system('reduce -BUILD -FLIP -Quiet {} > {}'.format(
                pdbpath, db_H_dir+vdm_dir+'/'+cluster+'/'+pdb))
