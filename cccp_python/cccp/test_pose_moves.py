import os
import numpy as np

import pyrosetta
from pyrosetta import rosetta

import loop_helix_loop_reshaping as LHLR

from loop_helix_loop_reshaping import simple_pose_moves
from loop_helix_loop_reshaping import pose_analysis


from generateCrickBB import *

#----------------------------
def get_rot_trans(mob_coords, targ_coords):
    mob_coords_com = mob_coords.mean(0)
    targ_coords_com = targ_coords.mean(0)
    mob_coords_cen = mob_coords - mob_coords_com
    targ_coords_cen = targ_coords - targ_coords_com
    cov_matrix = np.dot(mob_coords_cen.T, targ_coords_cen)
    U, S, Wt = np.linalg.svd(cov_matrix)
    R = np.dot(U, Wt)
    if np.linalg.det(R) < 0.:
        Wt[-1] *= -1
        R = np.dot(U, Wt)
    return R, mob_coords_com, targ_coords_com

#----------------------------
chains = 4
chL = 28
r0 = 7.36
r1 = 2.26
w0 = -2.45
w1 = 102.857
a = -12.01
ph1 = np.array([-9, -9, -9, -9])
cr = np.array([0, 1, 0])
zoff = np.array([0.0, 0.0, 0.0])
dph0 = np.array([90, 180, 270])
varargin = 'registerzoff'

XYZ = generateCrickBB(chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin)
xyz = XYZ[0:28]
#---------------------------

pyrosetta.init()

pose = pyrosetta.rosetta.core.pose.Pose()

pyrosetta.rosetta.core.pose.make_pose_from_sequence(pose, "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "fa_standard")

pose.annotated_sequence()

pose.size()

#Is it better to calculate the averange phi, psi and omega.
for i in range(1, 29):
    pose.set_phi(i, -64)
    pose.set_psi(i, -40)
    pose.set_omega(i, 180)

#output_dir = '/mnt/e/GitHub_Design/cccp_python'
#pose.dump_pdb(os.path.join(output_dir, 'test_unfinished.pdb'))
std_xyz = np.zeros((28, 3), dtype = float)
for i in range(1, 29):
    std_xyz[i-1] = [x for x in pose.residue(i).xyz(2)]

R, m_com, t_com = get_rot_trans(xyz, std_xyz)

xyz_transformed = np.dot((xyz - m_com), R) + t_com


xyz_v = pyrosetta.rosetta.numeric.xyzVector_double_t()
shift_xyz = pyrosetta.rosetta.numeric.xyzVector_double_t()
atomid = pyrosetta.rosetta.core.id.AtomID()

for i in range(1, 29):
    xyz = xyz_transformed[i-1,:]
    xyz_v.assign(xyz[0], xyz[1], xyz[2])
    atomid.set(2, i)
    shift_xyz.assign(xyz_v - pose.xyz(atomid))
    pose.set_xyz(atomid, xyz_v)
    for j in range(1, pose.residue(i).natoms() + 1):
        if j ==2:
            continue
        atomid.set(j, i)
        pose.set_xyz(atomid, pose.xyz(atomid) + shift_xyz)


output_dir = '/mnt/e/GitHub_Design/cccp_python'
pose.dump_pdb(os.path.join(output_dir, 'test_rt.pdb'))

#-----------------------------

sfxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function('score3')
sfxn.set_weight(rosetta.core.scoring.omega, 1)
sfxn.set_weight(rosetta.core.scoring.rama_prepro, 1)

mm = rosetta.core.kinematics.MoveMap()
for i in range(1, 29):
    mm.set_bb(i, True)

min_opts = rosetta.core.optimization.MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, True )
min_mover = rosetta.protocols.minimization_packing.MinMover()
min_mover.movemap(mm)
min_mover.min_options(min_opts)

min_mover.score_function(sfxn)
min_mover.apply(pose)

pose.dump_pdb(os.path.join(output_dir, 'test_rt_min.pdb'))
#-----------------------------
#Calculate phi psi, no rules find.
import prody as pr

import sys
sys.path.append(r'/mnt/e/GitHub_Design/smallprot')
from smallprot.struct_analysis import *


cccp_4_ala = '/mnt/e/GitHub_Design/cccp_python/std4H.pdb'
phi_180, psi_180, seq = meature_phipsi(cccp_4_ala)

protein = pr.parsePDB(cccp_4_ala)
om = []
for p in protein.iterResidues():
    try:
        om.append(pr.calcOmega(p, dist=None))
    except:
        om.append(None)

#-----------------------------
# Try to creat CaToAllAtom mover doesn't working.
xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
'''
<MOVERS>
    <CaToAllAtom name="ca2all_mover" />
</MOVERS>
''')
ca2all_mover = xmlobj.get_mover('ca2all_mover')


'''
#Example of creating MOVERS works.
fast_design_repeats=1
xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
'''
<MOVERS>
    <FastDesign name="fastdes" clear_designable_residues="0" repeats="{0}" ramp_down_constraints="1"/>
    <RotamerTrialsMover name="rot_trial" />
</MOVERS>
'''.format(fast_design_repeats))
fast_design = xmlobj.get_mover('fastdes')
rot_trial = xmlobj.get_mover('rot_trial')
'''
#------------------------------

pose = rosetta.core.import_pose.pose_from_file('std4H.pdb')

atomid.set(2, 2) 
x = pose.xyz(atomid) 
atomid.set(3, 2) 
y = pose.xyz(atomid)
print(x)
print(y)
print((x-y)[2])

#----------------------------


