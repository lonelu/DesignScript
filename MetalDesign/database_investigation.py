import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
import cluster
import transformation
import numpy as np

@dataclass
class pdb_info:
    pdb: str
    Metal: str
    validation: bool 
    pml_family:str
    clust_num: int
    selection: int

    def to_tab_string(self):
        pdb_info_str = self.pdb + '\t' + self.Metal + '\t' + str(self.validation) + '\t' + self.pml_family + '\t'+ str(self.clust_num) + '\t'+ str(self.selection)
        return pdb_info_str

# Handle metal

def check_metal_binding_validation(file_path, metal_sel):
    '''
    return if contain certain metal, if the metal is valid
    '''
    ### print('This is a pdb file')
    try:
        pdb_prody = pr.parsePDB(file_path)
    except:
        print('not sure')
        return False, False

    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return False, False
    for ni in nis:
        ni_index = ni.getIndex()
        all_near = pdb_prody.select('not water and within 2.6 of index ' + str(ni_index))
        if len(all_near)-1 >= 3:
            return True, True
    return True, False

# Write Summary

def write_pdb_info(filename, all_pdb_infos):
    '''
    Write information of all pdb.
    @ all_pdb_infos: [pdb_info]
    '''
    with open(filename, 'w') as f:
        f.write('pdb\tmetal\tvalidation\tpml_family\tclust_num\tselection\n')
        for r in all_pdb_infos:
            f.write(r.to_tab_string() + '\n')  

# Manipulate pml files

def extract_pml_pdb_info(workdir, metal_sel):
    '''
    extract redundancy infomation from protCID provided .pml file.
    '''
    all_pdb_infos = []

    for pml_path in os.listdir(workdir):
        if not pml_path.endswith("_byDomain.pml"):
            continue

        ### pml_path = 'ABC_tran_NI_byDomain.pml'
        the_set = []
        with open(workdir + pml_path, 'r') as f:
            while True:
                line = f.readline()
                if not line:
                    break
                if 'load' in line:
                    the_set.append(line.split(',')[0][5:])

            select = 0
            for pdb in the_set:
                contain_metal, is_valid = check_metal_binding_validation(workdir + pdb, metal_sel)
                if is_valid:
                    select += 1                
                metal = metal_sel[-2:] if contain_metal else 'none'
                selection = select if is_valid else 0
                the_pdb_info = pdb_info(pdb= pdb, Metal= metal, validation= is_valid, pml_family = pml_path, clust_num = len(the_set), selection = selection)
                all_pdb_infos.append(the_pdb_info)
    return all_pdb_infos
                
def extract_all_valid_pdbs():
    '''
    The function is Not used anywhere.
    Extract all valid pfam file. Return the filepath list.
    '''
    valid_metal_bind_pdbs = []

    for pdb_path in os.listdir(path):
        if not pdb_path.endswith(".pfam"):
            continue
            if check_metal_binding_validation(path + pdb_path):
                ### print('This is highly likely to be a real metal binding protein!')
                valid_metal_bind_pdbs.append(pdb_path)
    return valid_metal_bind_pdbs

def compare_rcsb_file(workdir, all_pdb_infos, all_rcsb_file = 'all_rcsb.txt', write_rcsb_protcid = False):
    '''
    read the rcsb file, write the coresponding all_pdb_infos from protcid.
    write all_rcsb_protcid.txt file.
    return the ones not in all_pdb_infos
    '''
    pfam_dict = dict()
    for i in range(len(all_pdb_infos)):
        pfam_dict['"' + all_pdb_infos[i].pdb[0:4] + '"'] = i

    all_pdbs = []
    pdb_names = []
    pdb_not_in_protcid = []
    with open(workdir + all_rcsb_file, 'r') as f:
        for r in f.readlines():
            pdb_name = str.lower(r.split('\t')[0])
            #print(pdb_name)
            pdb_names.append(pdb_name)
            if pdb_name !='':
                if pdb_name in pfam_dict.keys():
                    ind = pfam_dict[pdb_name]
                    r = r[:-1] + '\t' +  all_pdb_infos[ind].to_tab_string() + '\n'
                    #print(r)
                else:
                    pdb_not_in_protcid.append(pdb_name)
                all_pdbs.append(r)
    if write_rcsb_protcid:
        with open(workdir + 'all_rcsb_protcid.txt', 'w') as f:
            for r in all_pdbs:
                f.write(r)
    return pdb_names, pdb_not_in_protcid

def compare_protcid_file(workdir, all_pdb_infos, all_rcsb_file = 'all_rcsb.txt'):
    '''
    write protcid_rcsb_info.txt file.
    '''
    rcsb_dict = dict()
    with open(workdir + all_rcsb_file, 'r') as f:
        for r in f.readlines():
            pdb_name = str.lower(r.split('\t')[0])
            #print(pdb_name)
            if pdb_name !='':
                rcsb_dict[pdb_name] = r
    with open(workdir + 'protcid_rcsb_info.txt', 'w') as f:
        for r in all_pdb_infos:
            key = '"' + r.pdb[0:4] + '"'
            if key in rcsb_dict.keys():
                r = r.to_tab_string() + '\t' + rcsb_dict[key]
            f.write(r)

def pfam2pdb(workdir, file, sub_workdir = ''):
    if not os.path.exists(workdir + sub_workdir):
        os.mkdir(workdir + sub_workdir)
    #TO DO: do we need to remove the sequence file.
    # a_file = open("sample.txt", "r")
    # lines = a_file. readlines()
    # a_file. close()

    #Just copy the file.
    shutil.copy(workdir + file, workdir + sub_workdir + os.path.basename(file) + '.pdb')

def copy_pfam_subfolder(workdir, sub_workdir = 'pfam2pdbs/', info_select = []):
    '''
    Copy all pfam file into .pdb file in a subfolder.
    '''
    if len(info_select) > 0:
        _pdb = [p.pdb for p in info_select] 

    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pfam"):
            continue
        if len(info_select) > 0 and pdb_path not in _pdb:
            continue
        ld.pfam2pdb(workdir, pdb_path, sub_workdir)

# Manipulate rcsb file

def extract_rcsb_pdb_info(workdir_rcsb, metal_sel, sel_pdb = False, pdb_in=[]):
    '''
    extract validation infomation from rcsb download pdb files.
    '''
    p_dict = dict()
    for i in range(len(pdb_in)):
        p_dict[pdb_in[i].split('"')[1]] = i

    pdb_infos = []
    _pdbs = []
    for pdb in os.listdir(workdir_rcsb):
        if not pdb.endswith(".pdb"):
            continue
        #print(pdb[:-4])
        if sel_pdb and pdb[:-4] not in p_dict:
            continue
        _pdbs.append(pdb)

        contain_metal, is_valid = check_metal_binding_validation(workdir_rcsb + pdb, metal_sel)         
        metal = metal_sel[-2:] if contain_metal else 'none'
        selection = 1 if is_valid else 0
        the_pdb_info = pdb_info(pdb= pdb, Metal= metal, validation= is_valid, pml_family = '', clust_num = 0, selection = selection)
        pdb_infos.append(the_pdb_info)

    return pdb_infos


