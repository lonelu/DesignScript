# Sample dpp structures with all angles and distance on His.
# Superimpose them on the 6 helix backbone and filter by position.
# COMBS search potential interacting vdms.

import os
import prody as pr
import numpy as np
from scipy.spatial.transform import Rotation
from sklearn.neighbors import NearestNeighbors
import math
import pickle

def ligand_rot_is_clash(lig, interMolClashSets, interclash_dist = 3.0):
    '''
    The ligand it self could clash after rotation.
    The idea is to check the manually defined interMolClashSets 
    Such as interMolClashSets = [(['O1'], ['C1', 'C5']), (['O2'], ['C1', 'C5', 'C6'])], which means the rotation may induce clash between O1 and C1/C5.
    return True if clash
    '''
    for clashSet in interMolClashSets:
        atoms_rot = lig.select('serial ' + ' '.join(clashSet[0])).getCoords()
        rest_coords = lig.select('serial ' +  ' '.join(clashSet[1])).getCoords()

        nbrs = NearestNeighbors(radius= interclash_dist).fit(atoms_rot)
        adj_matrix = nbrs.radius_neighbors_graph(rest_coords).astype(bool)

        if np.sum(adj_matrix) >0:
            return True
    
    return False

def lig_dihedral_filter(lig, dh_atom_sel, dh_ranges):
    if len(dh_atom_sel)==0:
        return True
    dh = pr.calcDihedral(lig.select('serial ' + dh_atom_sel[0]), lig.select('serial ' + dh_atom_sel[1]), lig.select('serial ' + dh_atom_sel[2]), lig.select('serial ' + dh_atom_sel[3]))
    for r in dh_ranges:
        if dh > r[0] and dh < r[1]:
            return True
    return False

def rotate_ligs(orign_lig, rot, rest, rotation_degree, interMolClashSets, interclash_dist, dh_atom_sel, dh_ranges):

    lig = orign_lig.copy()

    # transform the lig to z axis.
    rot_coords = lig.select('serial ' + ' '.join(rot)).getCoords()
    rot_dist = pr.calcDistance(rot_coords[1], rot_coords[0])

    z_coords = np.zeros((2, 3), dtype=float)
    z_coords[1, -1] = rot_dist

    pr.calcTransformation(rot_coords, z_coords).apply(lig)

    all_ligs = []
    for i in range(0, 360, rotation_degree):
        _lig = lig.copy()
        rotation = Rotation.from_rotvec(np.radians(i)*np.array([0, 0, 1]))        
        _coords = rotation.apply(_lig.getCoords())
        _lig.setCoords(_coords)
        if len(rest) > 0:
            _lig.select('serial ' + ' '.join(rest)).setCoords(lig.select('serial ' + ' '.join(rest)).getCoords())

        pr.calcTransformation(_lig.select('serial ' + ' '.join(rest)).getCoords(), orign_lig.select('serial ' + ' '.join(rest)).getCoords()).apply(_lig)
        _lig.setTitle(lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i))
        #pr.writePDB(workdir + 'ligand_rotation/' +_lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i), _lig)
        if ligand_rot_is_clash(_lig, interMolClashSets, interclash_dist):
            continue
        if not lig_dihedral_filter(_lig, dh_atom_sel, dh_ranges):
            continue
        all_ligs.append(_lig)

    return all_ligs


def shift_lig(orign_lig, shift, rest, shift_bottom, shift_top, shift_step):

    lig = orign_lig.copy()

    # transform the lig to z axis.
    shift_coords = lig.select('serial ' + ' '.join(shift)).getCoords()
    shift_dist = pr.calcDistance(shift_coords[1], shift_coords[0])

    z_coords = np.zeros((2, 3), dtype=float)
    z_coords[1, -1] = shift_dist

    pr.calcTransformation(shift_coords, z_coords).apply(lig)

    all_ligs = []
    for i in range(shift_step + 1):
        direct = shift_bottom + (shift_top - shift_bottom)*i/(shift_step + 1) - shift_dist

        _lig = lig.copy()
        
        _coords = _lig.getCoords() 
        _coords[:, 2] += direct 
        _lig.setCoords(_coords)
        if len(rest) > 0:
            _lig.select('serial ' + ' '.join(rest)).setCoords(lig.select('serial ' + ' '.join(rest)).getCoords())

        pr.calcTransformation(_lig.select('serial ' + ' '.join(rest)).getCoords(), orign_lig.select('serial ' + ' '.join(rest)).getCoords()).apply(_lig)
        _lig.setTitle(lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i))
        #pr.writePDB(workdir + 'ligand_rotation/' +_lig.getTitle() + '_' + '-'.join(rot) + '_' + str(i), _lig)

        all_ligs.append(_lig)

    return all_ligs



workdir = '/mnt/e/DesignData/Metalloprotein/DPP/'
outdir = workdir + 'test/'

ligfile = workdir + 'fe-dpp-cooh.pdb'

orign_lig = pr.parsePDB(ligfile)

#lig.getSerials() 


rot = ['7', '69']
rest = ['2', '3', '4', '5', '6', '7', '8', '1', '9', '10', '11', '12']
dh_atom_sel = ['6', '7', '69', '21']
dh_range = [(-93.7, -33.7)]


rot1 = ['21', '37']
rest1_r = ['73', '61', '75', '74', '40', '41', '62', '42', '63', '37', '38', '59', '39', '60']
rest1 = [str(x) for x in range(1, 76) if str(x) not in rest1_r]
interMolClashSets1 = [(['50', '51'], ['59', '38', '42', '63'])]
interclash_dist1 = 2.7

rot2 = ['31', '43']
rest2_r = ['43', '48', '68', '47', '67', '46', '66', '70', '72', '71', '45', '65', '44', '64']
rest2 = [str(x) for x in range(1, 76) if str(x) not in rest2_r]
interMolClashSets2 = [(['55', '56'], ['48', '68', '44', '64'])]
interclash_dist2 = 2.7

shift = ['7', '69']
shift_bottom = 2.0
shift_top = 3.0
shift_step = 4



def generate_rotated_ligs():
    '''
    The method only works for 2 rots.
    TO DO: use recursive algorithm to allow more than 2 rots.
    '''
    all_ligs = []
    i = 0
    ligs = rotate_ligs(orign_lig, rot, rest, 2, [], 0, dh_atom_sel, dh_range)
    for _lig in ligs:
        ligs1 = rotate_ligs(_lig, rot1, rest1, 18, interMolClashSets1, interclash_dist1, [], [])
        for _lig1 in ligs1:
            ligs2 = rotate_ligs(_lig1, rot2, rest2, 18, interMolClashSets2, interclash_dist2, [], [])
            for _lig2 in ligs2:
                ligs3 = shift_lig(_lig2, shift, rest, shift_bottom, shift_top, shift_step)
                all_ligs.extend(ligs3)
    return all_ligs


all_ligs = generate_rotated_ligs()
len(all_ligs)

#for i in range(100):
#    lg = all_ligs[i*5]
#    pr.writePDB(outdir + lg.getTitle() + '.pdb', lg)

with open(outdir + 'all_ligs.pkl', 'wb') as f:
    pickle.dump(all_ligs, f)

with open(outdir + 'all_ligs.pkl', 'rb') as f:
    alligs = pickle.load(f)