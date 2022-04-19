# This python script take a helix bundle pose, and out put a pdb with defined TM region with PDB label info and a .span file
# This python script is written by Hua Bai and Peilong Lu.
# This python is python2
# huabai@uw.edu lpl15@uw.edu

import math
import pyrosetta
import argparse
# from datetime import datetime
from glob import glob


## Step1: directly check per residue sasa score
def based_on_sasa(input_pose, ball_size=1, sasa_threshold = 10):
    ## define two vector placeholders for the scores. One for atoms, one for residues
    atom_sasa = pyrosetta.rosetta.core.id.AtomID_Map_double_t()
    rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
    pyrosetta.rosetta.core.scoring.calc_per_atom_sasa(input_pose,atom_sasa,rsd_sasa,ball_size)
    selected_by_sasa = []
    for idx, resi_sasa in enumerate(rsd_sasa):
        if resi_sasa > sasa_threshold:
            selected_by_sasa.append(idx+1)
    return selected_by_sasa


## Step2: select polar residues
def based_on_polarity(input_pose):
    selected_by_polar = []
    for i in range(input_pose.size()):
        idx = i+1
    if input_pose.residue_type(idx).is_polar():
        selected_by_polar.append(idx)
    return selected_by_polar


## Step2.1: select chain A
def chain_A(input_pose):
    selected_by_chain = []
    chaina = input_pose.split_by_chain()[1]
    for i in range(chaina.size()):
        idx = i+1
        selected_by_chain.append(idx)
    return selected_by_chain


## Step3-Method1: select outward pointing residues by angles. NOT USE.
## output unit: degree
## the angle between CA-CB ray and origin-CA ray
## for example, if the angle larger than 90 degree, smaller than 270 degree, the CA-CB may point to inward
def calc_angle(pose, resi_idx):
    ## special case: GLY
    if pose.residue(resi_idx).name() == "GLY":
        return 999

    ## for CA
    ca_x = pose.residue(resi_idx).atom("CA").xyz()[0]
    ca_y = pose.residue(resi_idx).atom("CA").xyz()[1]
    ca_angle = math.atan2(ca_y,ca_x)
    ## for CB
    cb_x = pose.residue(resi_idx).atom("CB").xyz()[0]
    cb_y = pose.residue(resi_idx).atom("CB").xyz()[1]
    cb_a_x = cb_x - ca_x
    cb_a_y = cb_y - ca_y
    cb_a_angle = math.atan2(cb_a_y,cb_a_x)
    ## return the difference
    return ((cb_a_angle - ca_angle) / math.pi * 180) % 360


def based_on_CA_CB_angle(input_pose, angle_threshold = 120):
    selected_by_angle = []
    for resi_idx in range(input_pose.size()):
        check_angle = calc_angle(input_pose, resi_idx+1)
        if check_angle < angle_threshold or check_angle > 360-angle_threshold and check_angle != 999:
            selected_by_angle.append(resi_idx+1)
    return selected_by_angle


## Step3-Method2: select outward pointing residues by distance differencs
## output unit: distance, A
## the compare the distances between origin-CA and both origin-N ray and origin-C ray
## sum the distance differences
## the more negative the difference, the higher chance it is pointing inward
def compare_distance_to_origin(pose, resi_idx):
    n_x = pose.residue(resi_idx).atom("N").xyz()[0]
    n_y = pose.residue(resi_idx).atom("N").xyz()[1]
    # n_z = pose.residue(resi_idx).atom("N").xyz()[2]
    n_distance = math.sqrt(n_x**2 + n_y**2)
    c_x = pose.residue(resi_idx).atom("C").xyz()[0]
    c_y = pose.residue(resi_idx).atom("C").xyz()[1]
    # c_z = pose.residue(resi_idx).atom("C").xyz()[2]
    c_distance = math.sqrt(c_x**2 + c_y**2)
    ca_x = pose.residue(resi_idx).atom("CA").xyz()[0]
    ca_y = pose.residue(resi_idx).atom("CA").xyz()[1]
    # ca_z = pose.residue(resi_idx).atom("CA").xyz()[2]
    ca_distance = math.sqrt(ca_x**2 + ca_y**2)
    ## CA-N distance, corrected by helix pitch angle, not use
    ca_n_distance_diff = ca_distance - n_distance
    ## ca_n_z_diff = ca_z - n_z
    ## angle = math.atan2(ca_n_z_diff, math.fabs(ca_n_distance_diff))
    ## ca_n_distance_diff_corrected = ca_n_distance_diff / math.cos(angle)
    ## CA-C distance, corrected by helix pitch angle, not use
    ca_c_distance_diff = ca_distance - c_distance
    ## ca_c_z_diff = ca_z - c_z
    ## angle = math.atan2(ca_c_z_diff, math.fabs(ca_c_distance_diff))
    ## ca_c_distance_diff_corrected = ca_c_distance_diff / math.cos(angle)
    return ca_n_distance_diff + ca_c_distance_diff


def based_on_distance(input_pose, dist_diff_threshold = -0.6):
    selected_by_distance = []
    for resi_idx in range(input_pose.size()):
        check_distance = compare_distance_to_origin(input_pose, resi_idx+1)
        if check_distance > dist_diff_threshold:
            selected_by_distance.append(resi_idx+1)
    return selected_by_distance


## Step4: choose residues inside the membrane
import numpy as np
def survey_z(input_pose, percentile , midpoint):
    all_z = []
    for resi_idx in range(input_pose.size()):
        ca_z = input_pose.residue(resi_idx+1).atom("CA").xyz()[2]
        all_z.append(ca_z)
    N_z=all_z[1]
    C_z=all_z[-1]
    z_array = np.array(all_z)
    selected_by_z = []
    for idx, value in enumerate(z_array):
        if float(value) > float(midpoint-percentile) and float(value) < float(midpoint+percentile):
            selected_by_z.append(idx+1)
    return selected_by_z


def WY_z(input_pose, percentile , midpoint):
    all_z = []
    for resi_idx in range(input_pose.size()):
        ca_z = input_pose.residue(resi_idx+1).atom("CA").xyz()[2]
        all_z.append(ca_z)
    N_z=all_z[1]
    C_z=all_z[-1]
    z_array = np.array(all_z)
    selected_by_z = []
    for idx, value in enumerate(z_array):
        if float(value) > float(midpoint+percentile-2.5) and float(value) < float(midpoint+percentile+2.5):
            selected_by_z.append(idx+1)
    return selected_by_z


def RK_z(input_pose, percentile , midpoint):
    all_z = []
    for resi_idx in range(input_pose.size()):
        ca_z = input_pose.residue(resi_idx+1).atom("CA").xyz()[2]
        all_z.append(ca_z)
    N_z=all_z[1]
    C_z=all_z[-1]
    z_array = np.array(all_z)
    selected_by_z = []
    for idx, value in enumerate(z_array):
        if float(value) > float(midpoint-percentile-2.5) and float(value) < float(midpoint-percentile+2.5):
            selected_by_z.append(idx+1)
    return selected_by_z



def F_z(input_pose, percentile , midpoint):
    all_z = []
    for resi_idx in range(input_pose.size()):
        ca_z = input_pose.residue(resi_idx+1).atom("CA").xyz()[2]
        all_z.append(ca_z)
    N_z=all_z[1]
    C_z=all_z[-1]
    z_array = np.array(all_z)
    selected_by_z = []
    for idx, value in enumerate(z_array):
        if float(value) > float(midpoint-1.2) and float(value) < float(midpoint+1.2):
            selected_by_z.append(idx+1)
    return selected_by_z


"""import numpy as np
def survey_z(input_pose, percentile = 90):
    all_z = []
    for resi_idx in range(input_pose.size()):
        ca_z = input_pose.residue(resi_idx+1).atom("CA").xyz()[2]
        all_z.append(ca_z)
    z_array = np.array(all_z)
    selected_by_z = []
    for idx, value in enumerate(z_array):
        if value > np.percentile(z_array, (100-percentile)/2) and value < np.percentile(z_array, 100- (100-percentile)/2):
            selected_by_z.append(idx+1)
        return selected_by_z
"""

## Step5: find the interdections of the selected residues
## https://stackoverflow.com/questions/3852780/python-intersection-of-multiple-lists
def intersect(*d):
    sets = iter(map(set, d))
    result = sets.next()
    for s in sets:
        result = result.intersection(s)
    return result



## Step6: Combine the previous steps, and add remarks to pdb file
from shutil import copyfile
import os
## pdb_input: the location/path of the input pdb
## final_selection: the list of index of the selected residues, 1-based
def add_PDBinfo_remark (pdb_input, final_selection,WY_selection,RK_selection,F_selection):
    path_list = list(os.path.split(pdb_input))
    design_id = path_list[-1].split(".pdb")[0]
    modified_file_name = "{}_with_remarks.pdb".format(design_id)
    path_list[-1] = modified_file_name
    output_path = os.path.join(*path_list)
    try:
        os.remove(output_path)
    except OSError:
        pass
    os.system("cat "+pdb_input+"|~/scripts/filt_chain.pl A > "+output_path)
    with open(output_path,"a") as output_pdb:
        for idx in final_selection:
            pdbinfo_label_remark = "REMARK PDBinfo-LABEL: {:>4} OUTWARD_POLAR\n".format(idx)
            output_pdb.write(pdbinfo_label_remark)
        for idx in WY_selection:
            pdbinfo_label_remark = "REMARK PDBinfo-LABEL: {:>4} WY_RING\n".format(idx)
            output_pdb.write(pdbinfo_label_remark)
        for idx in RK_selection:
            pdbinfo_label_remark = "REMARK PDBinfo-LABEL: {:>4} RK_RING\n".format(idx)
            output_pdb.write(pdbinfo_label_remark)
        for idx in F_selection:
            pdbinfo_label_remark = "REMARK PDBinfo-LABEL: {:>4} F_RING\n".format(idx)
            output_pdb.write(pdbinfo_label_remark)


"""
## Step6.1: add to run.sh
from shutil import copyfile
import os
## pdb_input: the location/path of the input pdb
## final_selection: the list of index of the selected residues, 1-based
def add_PDBinfo_remark (pdb_input, final_selection,WY_selection,RK_selection):
    path_list = list(path.split(pdb_input))
    design_id = path_list[-1]
    modified_file_name = "{}.sh".format(design_id)
    path_list[-1] = modified_file_name
    output_path = os.path.join(*path_list)
    try:
        os.remove(output_path)
    except OSError:
        pass
    os.system("cat "+pdb_input+"|~/scripts/filt_chain.pl A > "+output_path)
    with open(output_path,"a") as output_pdb:
        for idx in final_selection:
            pdbinfo_label_remark = "REMARK PDBinfo-LABEL: {:>4} OUTWARD_POLAR\n".format(idx)
            output_pdb.write(pdbinfo_label_remark)
        for idx in WY_selection:
            pdbinfo_label_remark = "REMARK PDBinfo-LABEL: {:>4} WY_RING\n".format(idx)
            output_pdb.write(pdbinfo_label_remark)
        for idx in RK_selection:
            pdbinfo_label_remark = "REMARK PDBinfo-LABEL: {:>4} RK_RING\n".format(idx)
            output_pdb.write(pdbinfo_label_remark)
"""


## create span file for illustration in pymol
def add_span_file(selection_keyword, selection_list, length):
    with open("{}.span".format(selection_keyword), "w") as script:
        script.write("TM region prediction for BRD4 predicted using OCTOPUS\n")
        TM_split=[]
        for idx in range(len(selection_list)-1):
            if selection_list[idx] != selection_list[idx+1]-1:
                TM_split.append(idx)
        TM_N = len(TM_split)+1
        script.write("{} {}\n".format(TM_N, length))
        script.write("antiparallel\n")
        script.write("n2c\n")
        for i in range(TM_N):
            if i ==0:
                script.write("{} {} {} {}\n".format(selection_list[0],selection_list[TM_split[i]],selection_list[0],selection_list[TM_split[i]] ))
            elif i == TM_N-1:
                script.write("{} {} {} {}\n".format(selection_list[TM_split[-1]+1],selection_list[-1],selection_list[TM_split[-1]+1],selection_list[-1] ))
            else:
                script.write("{} {} {} {}\n".format(selection_list[TM_split[i-1]+1],selection_list[TM_split[i]],selection_list[TM_split[i-1]+1],selection_list[TM_split[i]] ))


## create pml file for illustration in pymol
def create_pml_file(selection_keyword, selection_list):
    with open("select_surface_{}.pml".format(selection_keyword), "w") as script:
        script.write("select {}, resi ".format(selection_keyword))
    for idx, resi_idx in enumerate(selection_list):
        if idx == 0:
            script.write("{}".format(resi_idx))
        else:
            script.write("+{}".format(resi_idx))
    script.write("\n")
    script.write("util.cbag all\n")
    script.write("util.cbay ({})\n".format(selection_keyword))


## argparse argument
def get_argparse():
    parser = argparse.ArgumentParser(description='Select the surface polar residues for memebrane channel designs')
    parser.add_argument("-b", '--ball_size', metavar = "\b", type=float, dest='ball_size',
    default=1,
    help='The ball size for calculating sasa, the smaller the ball, the deeper the surface can reach. The default value is 1')
    parser.add_argument("-s", '--sasa_threshold', metavar = "\b", type=float, dest='sasa_threshold',
    default=10,
    help='The sasa threshold for determining whether the residue is surface or core. The higher the sasa, the residue is more like to be on surface. The default value is 10')
    parser.add_argument("-d", '--distance_threshold', metavar = "\b", type=float, dest='distance_threshold',
    default=-0.6,
    help='the sum of the distance differences of (CA-axis - N-axis) and (CA-axis - C-axis), the more negative the sum is, the more likely the residue is pointing inward. The default value is -0.6')
    parser.add_argument("-z", '--z_percentile', metavar = "\b", type=float, dest='z_percentile',
    default=13.5,
    help='the percentage of the residues are considered to be inside membrane. The default value is 3')
    parser.add_argument("-i", '--input_dir', metavar = "\b", type=str, dest='input_dir',
    default=".",
    help='This is where the input pdbs are located. The default is current dir "."')
    parser.add_argument("-m", '--imidpoint', metavar = "\b", type=float, dest='imidpoint',
    default=0,
    help='the midpoint in the membrane. The default value is 0')
    ## parser.add_argument('--output_dir', type=str, dest='output_dir',
    ## default=datetime.now().strftime("pdbs_with_remarks_%y%m%d_%p_%I%M%S"),
    ## help='This is where the pdbs with remarks will be stored')
    return parser

if __name__ == "__main__":
    ## load arguments
    parser = get_argparse()
    args=parser.parse_args("-m -10 -z 15 -s 8".split()) #This is important for running in jupyter
    # args=parser.parse_args()
    ## init pyrosetta
    pyrosetta.init(options="")
    target_pdbs = glob("{}/*4.pdb".format(args.input_dir))
    for target in target_pdbs:
        path_list = os.path.split(target)
        design_id = path_list[-1].split(".pdb")[0]
        input_pose = pyrosetta.pose_from_file(target)
        length = len(chain_A(input_pose))
        final_selection = list(intersect(based_on_sasa(input_pose, args.ball_size, args.sasa_threshold), based_on_polarity(input_pose),
        based_on_distance(input_pose, args.distance_threshold), chain_A(input_pose), survey_z(input_pose, args.z_percentile,
        args.imidpoint) ))
        WY_selection = list(intersect(based_on_sasa(input_pose, args.ball_size, args.sasa_threshold), based_on_polarity(input_pose),
        based_on_distance(input_pose, args.distance_threshold), chain_A(input_pose), WY_z(input_pose, args.z_percentile,
        args.imidpoint) ))
        RK_selection = list(intersect(based_on_sasa(input_pose, args.ball_size, args.sasa_threshold), based_on_polarity(input_pose),
        based_on_distance(input_pose, args.distance_threshold), chain_A(input_pose), RK_z(input_pose, args.z_percentile,
        args.imidpoint) ))
        F_selection = list(intersect(based_on_sasa(input_pose, args.ball_size, args.sasa_threshold), based_on_polarity(input_pose),
        based_on_distance(input_pose, args.distance_threshold), chain_A(input_pose), F_z(input_pose, args.z_percentile,
        args.imidpoint) ))
        TM_selection = list(intersect(chain_A(input_pose), survey_z(input_pose, args.z_percentile, args.imidpoint)))
        add_span_file("{}".format(design_id), TM_selection, length)
        add_PDBinfo_remark(target,final_selection,WY_selection,RK_selection,F_selection)
        create_pml_file("{}_final".format(design_id), final_selection)
        create_pml_file("{}_RK".format(design_id), WY_selection)
        create_pml_file("{}_WY".format(design_id), RK_selection)
        create_pml_file("{}_F".format(design_id), F_selection)  