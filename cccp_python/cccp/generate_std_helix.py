import os
import pyrosetta
from pyrosetta import rosetta

pyrosetta.init()

pose = pyrosetta.rosetta.core.pose.Pose()

cc_len = 28
seq_ala = ''.join(['A' for a in list(range(cc_len))])

pyrosetta.rosetta.core.pose.make_pose_from_sequence(pose, seq_ala, "fa_standard")

for i in range(1, 29):
    pose.set_phi(i, -57)
    pose.set_psi(i, -47)
    pose.set_omega(i, 180)

output_dir = '/mnt/e/GitHub_Design/cccp_python'
pose.dump_pdb(os.path.join(output_dir, 'test_std_helix.pdb'))