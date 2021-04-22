# reload module in ipython
import sys
import importlib 

importlib.reload(sys.modules['ligand_database']) 
from ligand_database import * 

# rename files
import os

def change_name(workdir, metal):
    for file in os.listdir(workdir):
        if file.endswith(".pdb"):
            src = workdir + file
            dst = workdir + file.split('_')[0] + metal + file.split('_')[1]
            os.rename(src, dst)

workdir = '/mnt/e/DesignData/ligands/ZN_rcsb/1_seq_cores/'

metal = '_ZN_'

change_name(workdir, metal)

# copy files
import os 
import shutil

def get_pdbs(workdir):
    pdbs = []
    for file in os.listdir(workdir):
        if file.endswith(".pdb"):
            pdb = file.split('_')[4] + '_' + file.split('_')[5] + '_' + file.split('_')[6]
            pdbs.append(pdb)
    return pdbs

def copy_pdbs(workdir, outdir, pdbs):
    for file in os.listdir(workdir):
        if file.endswith(".pdb"):
            if file.split('.')[0] in pdbs:
                shutil.copy(workdir + file, outdir + file)

workdir = '/mnt/e/DesignData/ligands/all_metal/13_2his_cores_cluster_7/0/'

pdbs = get_pdbs(workdir)

workdir = '/mnt/e/DesignData/ligands/all_metal/2_seq_cores_reps/'

outdir = '/mnt/e/DesignData/ligands/all_metal/13_2his_cores_cluster_7_0_cores/'

copy_pdbs(workdir, outdir, pdbs)