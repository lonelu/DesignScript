import os
import pyrosetta
from pyrosetta import rosetta
import pyrosetta.toolbox.mutants as mutate

'''
python /mnt/e/GitHub_Design/DesignScript/pyrosetta/pyrosetta_mutation.py
'''
pyrosetta.init(extra_options="-ignore_zero_occupancy false ") 

def mutate_all_to_ala_or_gly(indir, outdir, pdb_name, aa = 'G'):
    os.makedirs(outdir, exist_ok=True)
    pose = rosetta.core.import_pose.pose_from_file(indir + pdb_name)
    for i in range(1, pose.pdb_info().nres()+ 1):
        mutate.mutate_residue(pose, i, aa)
    pose.dump_pdb(os.path.join(outdir, pdb_name[:-4] + '_' + aa + '.pdb'))
    return

def main():
    indir = '/mnt/e/DesignData/Metalloenzyme/HelixZn/'
    outdir = indir + 'all_gly/'

    for file in os.listdir(indir):
        if not '.pdb' in file:
            continue
        mutate_all_to_ala_or_gly(indir, outdir, file, aa = 'G')
    return

main()
