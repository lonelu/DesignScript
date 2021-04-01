import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil

@dataclass
class pdb_info:
    pdb: str
    Metal: str
    validation: bool 
    pml_family:str
    clust_num: int
    selection: int

def to_tab_string(pdb_info):
    pdb_info_str = pdb_info.pdb + '\t' + pdb_info.Metal + '\t' + str(pdb_info.validation) + '\t' + pdb_info.pml_family + '\t'+ str(pdb_info.clust_num) + '\t'+ str(pdb_info.selection) + '\n'
    return pdb_info_str

def check_metal_binding_validation(file_path, metal_sel):
    '''
    return if contain certain metal, if the metal is valid
    '''
    ### print('This is a pdb file')
    pdb_prody = pr.parsePDB(file_path)
    coords = pdb_prody.getCoords()
    ### coords = np.array(pdb_prody.getCoords(), dtype='float32') 

    # A pdb can contain more than one NI.
    if not pdb_prody.select(metal_sel):
        print(file_path)
        return False, False
    for ni in pdb_prody.select(metal_sel):
        ni_coords = ni.getCoords()
        ni_index = ni.getIndex()

        dists = cdist(coords, coords)
        ### dists[ni_index, ni_index] == 0

        # The rule for cutting is: the metal has equal or more than 3 'N' or 'O' atom within 2.6 A, 
        dist_lim_index =  [idx for idx, val in enumerate(dists[ni_index]) if val < 2.6 and val >0]

        all_near = []
        for dli in dist_lim_index:
            a = pdb_prody[dli]
            if a.getElement() == 'N' or a.getElement()=='O':
                all_near.append(a.getResname())

        if len(all_near) >= 3:
            return True, True
    return True, False
    
def pfam2pdb(workdir, file, sub_workdir):
    if not os.path.exists(workdir + sub_workdir):
        os.mkdir(workdir + sub_workdir)
    #TO DO: do we need to remove the sequence file.
    # a_file = open("sample.txt", "r")
    # lines = a_file. readlines()
    # a_file. close()

    #Just copy the file.
    shutil.copy(workdir + file, workdir + sub_workdir + os.path.basename(file) + '.pdb')


# extract redundancy infomation from protCID provided .pml file.
def extract_pml_pdb_info(workdir, metal_sel):
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
                
def write_pdb_info(filename, all_pdb_infos):
    '''
    Write information of all pdb.
    @ loop_infos: [pdb_info]
    '''
    with open(filename, 'w') as f:
        f.write('pdb\tmetal\tvalidation\tpml_family\tclust_num\tselection\n')
        for r in all_pdb_infos:
            f.write(r.pdb + '\t' + r.Metal + '\t' + str(r.validation) + '\t'
                + r.pml_family + '\t'+ str(r.clust_num) + '\t'+ str(r.selection) + '\n')  

def extract_all_valid_pdbs():
    '''
    The function is Not used anywhere.
    '''
    valid_metal_bind_pdbs = []

    for pdb_path in os.listdir(path):
        if not pdb_path.endswith(".pfam"):
            continue
            if check_metal_binding_validation(path + pdb_path):
                ### print('This is highly likely to be a real metal binding protein!')
                valid_metal_bind_pdbs.append(pdb_path)
    return valid_metal_bind_pdbs

def organize_rcsb_file(workdir = "/mnt/e/DesignData/ligands/NI_rcsb/"):
    all_lines = []
    for file in os.listdir(workdir):
        if file.endswith(".csv"):
            with open(workdir + file, 'r') as f:
                all_lines = f.readlines()
    with open(workdir + 'all_rcsb.txt', 'w') as f:
        f.write('\t'.join(all_lines[0].split(',')))
        for r in all_lines:
            if 'Entry ID' not in r and r.split(',')[0]!= '':
                f.write('\t'.join(r.split(',')))

def compare_rcsb_file(workdir, all_pdb_infos, all_rcsb_file = 'all_rcsb.txt'):
    '''
    read the rcsb file, write the coresponding all_pdb_infos from protcid.

    '''
    pfam_dict = dict()
    for i in range(len(all_pdb_infos)):
        pfam_dict['"' + all_pdb_infos[i].pdb[0:4] + '"'] = i

    all_pdbs = []
    pdb_names = []
    with open(workdir + all_rcsb_file, 'r') as f:
        for r in f.readlines():
            pdb_name = str.lower(r.split('\t')[0])
            #print(pdb_name)
            pdb_names.append(pdb_name)
            if pdb_name !='':
                if pdb_name in pfam_dict.keys():
                    ind = pfam_dict[pdb_name]
                    r = r[:-1] + '\t' +  to_tab_string(all_pdb_infos[ind])
                    print(r)
                all_pdbs.append(r)
    with open(workdir + 'all_rcsb_protcid.txt', 'w') as f:
        for r in all_pdbs:
            f.write(r)

def compare_protcid_file(workdir, all_pdb_infos, all_rcsb_file = 'all_rcsb.txt'):
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
                r = to_tab_string(r)[:-1] + '\t' + rcsb_dict[key]
            f.write(r)

workdir = "/mnt/e/DesignData/ligands/NI/"
metal_sel = 'name NI'

all_pdb_infos = extract_pml_pdb_info(workdir, metal_sel)

write_pdb_info(workdir + 'protcid_info.txt', all_pdb_infos)

# # copy all pfam file into .pdb file in a subfolder.
# for pdb_path in os.listdir(workdir):
#     if not pdb_path.endswith(".pfam"):
#         continue
#     pfam2pdb(workdir, pdb_path, 'pfam2pdbs/')
    
   








