import os
import numpy as np

import pyrosetta
from pyrosetta import rosetta

from _generateCrickBB import *

pyrosetta.init()

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
def generate_default_CrickBB()
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

    return xyz


def generate_cc(output_dir, xyz, cc_len = 28):

    seq_ala = ''.join(['A' for a in list(range(cc_len))])

    pose = pyrosetta.rosetta.core.pose.Pose()

    pyrosetta.rosetta.core.pose.make_pose_from_sequence(pose, seq_ala, "fa_standard")

    pose.annotated_sequence()

    pose.size()

    #Is it better to calculate the averange phi, psi and omega.
    for i in range(1, cc_len + 1):
        pose.set_phi(i, -64)
        pose.set_psi(i, -40)
        pose.set_omega(i, 180)

    #output_dir = '/mnt/e/GitHub_Design/cccp_python'
    #pose.dump_pdb(os.path.join(output_dir, 'test_unfinished.pdb'))
    std_xyz = np.zeros((cc_len, 3), dtype = float)
    for i in range(1, cc_len + 1):
        std_xyz[i-1] = [x for x in pose.residue(i).xyz(2)]

    R, m_com, t_com = get_rot_trans(xyz, std_xyz)

    xyz_transformed = np.dot((xyz - m_com), R) + t_com


    xyz_v = pyrosetta.rosetta.numeric.xyzVector_double_t()
    shift_xyz = pyrosetta.rosetta.numeric.xyzVector_double_t()
    atomid = pyrosetta.rosetta.core.id.AtomID()

    for i in range(1, cc_len + 1):
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

    pose.dump_pdb(os.path.join(output_dir, 'output.pdb'))
