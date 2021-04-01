import os
import pyrosetta
from pyrosetta import rosetta

pyrosetta.init()

#--------------------------pose from sequence
pose = pyrosetta.rosetta.core.pose.Pose()
pyrosetta.rosetta.core.pose.make_pose_from_sequence(pose, "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "fa_standard")


#--------------------------pose from file
## If want to add ignore_zero_occupancy
## One can add the parameter '-ignore_zero_occupancy' in the init step. 
pyrosetta.init(extra_options="-ignore_zero_occupancy false ") 

input_pdb = '/mnt/e/GitHub_Design/scripts/data/sasa/full_7a_14a.pdb'
pose = rosetta.core.import_pose.pose_from_file(input_pdb)

## Or use a very verbose way. ## TO CHECK: maybe there is a better way.
opt = pyrosetta.rosetta.core.import_pose.ImportPoseOptions()
opt.set_ignore_zero_occupancy(False)
input_pdb_4h = '/mnt/e/GitHub_Design/scripts/data/sasa/full_7a_10a.pdb'
pose_4h = pyrosetta.rosetta.core.pose.Pose()
pyrosetta.rosetta.core.import_pose.pose_from_file(pose_4h, input_pdb_4h, opt)

#---------------------------save pose to pdb
import os
output_dir = '/mnt/e/GitHub_Design/cccp_python'
pose.dump_pdb(os.path.join(output_dir, 'test_std_helix.pdb'))

#---------------------------calculate Sasa.
sa = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
sa.calculate(pose)
#sa.calculate(pose_4h)
sa.get_atom_sasa()

## TO DO: the following is not working.
# sa_le = pyrosetta.rosetta.core.scoring.sasa.LeGrandSasa()
# sa_le.calculate(pose)
# sa_le.calculate(pose_4h)

