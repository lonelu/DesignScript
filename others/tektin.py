import os
from pyrosetta import *
from pyrosetta import rosetta
import pyrosetta.toolbox.mutants as mutate
pyrosetta.init(extra_options="-ignore_zero_occupancy false ") 

#--------------------------pose from sequence
input_pdb = '/mnt/e/DesignData/_temp/ZhenChen/tektin4_rb.pdb'
pose = rosetta.core.import_pose.pose_from_file(input_pdb)

from Bio import SeqIO
input_fasta = '/mnt/e/DesignData/_temp/ZhenChen/tektin5_rb.fasta'
fasta_sequence = list(SeqIO.parse(open(input_fasta),'fasta'))[0]


for i in range(1, len(fasta_sequence.seq)+ 1):
    if pose.sequence()[i-1] == fasta_sequence.seq[i-1]:
        continue
    mutate.mutate_residue(pose, i, fasta_sequence.seq[i-1])

output_dir = '/mnt/e/DesignData/_temp/ZhenChen'
pose.dump_pdb(os.path.join(output_dir, 'tektin_4_muta_5.pdb'))

scorefxn = rosetta.core.scoring.get_score_function()
mm = MoveMap()
mm.set_bb(False)
mm.set_chi(True)

def min_pose(pose, scorefxn, mm):
    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover.movemap(mm)
    min_mover.score_function(scorefxn)
    min_mover.min_type("lbfgs_armijo")
    min_mover.tolerance(1e-6)
    min_mover.apply(pose)

min_pose(pose,scorefxn,mm)

output_dir = '/mnt/e/DesignData/_temp/ZhenChen'
pose.dump_pdb(os.path.join(output_dir, 'tektin_4_muta_5_minmover.pdb'))

