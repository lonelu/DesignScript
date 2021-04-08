import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Qbits')
from qbits import cluster, transformation
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

def organize_rcsb_file(workdir = "/mnt/e/DesignData/ligands/NI_rcsb/"):
    '''
    The .csv files downloaded from rcsb database will be combined first, 
    then generate tab deliminated txt file.
    '''
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

# Extract peptide

def connectivity_filter(pdb_prody, ind, ext_ind):
    res1 = pdb_prody.select('resindex ' + str(ind))
    res2 = pdb_prody.select('resindex ' + str(ext_ind))
    if not res2:
        return False
    if res1[0].getResnum() - res2[0].getResnum() == ind - ext_ind and res1[0].getChid() == res2[0].getChid() and res1[0].getSegname() == res2[0].getSegname():
            return True
    return False

def extend_res_indices(inds_near_res, pdb_prody, extend = 2):
    extend_inds = []
    inds = set()
    for ind in inds_near_res:
        for i in range(-2, 3):
            the_ind = ind + i
            if the_ind not in inds and the_ind>= 0 and connectivity_filter(pdb_prody, ind, the_ind):
                extend_inds.append(the_ind)
                inds.add(the_ind)
    return extend_inds

def get_metal_core_seq(file_path, metal_sel):
    try:
        pdb_prody = pr.parsePDB(file_path)
    except:
        print('not sure')
        return []

    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return []
    
    metal_cores = []
    count = 0
    for ni in nis:
        ni_index = ni.getIndex()
        all_near = pdb_prody.select('nitrogen or oxygen').select('not water and within 2.6 of index ' + str(ni_index))
        if not all_near or len(all_near) < 3:
            continue          
        inds = all_near.getResindices()
        all_near_res = pdb_prody.select('protein and resindex ' + ' '.join([str(ind) for ind in inds]))
        if not all_near_res or len(all_near_res) < 1:
            continue     
        inds_near_res =  all_near_res.getResindices()
        ext_inds = extend_res_indices(inds_near_res, pdb_prody)
        count += 1
        sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
        metal_cores.append((os.path.basename(file_path).split('.')[0] + '_' + str(count), sel_pdb_prody))        
    return metal_cores

def extract_all_core_seq(workdir, metal_sel):
    cores = []
    first_get = False
    first = None

    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        core = get_metal_core_seq(workdir + pdb_path, metal_sel)
        if not first_get:
            first = core[0]
            first_get = True
        cores.extend(core)
    return cores, first

def superimpose_core_and_writepdb(cores, first, workdir, metal_sel, sub_workdir = 'seq_cores/'):
    if not os.path.exists(workdir + sub_workdir):
        os.mkdir(workdir + sub_workdir)

    for c in cores:
        a_coords = first[1].select(metal_sel)[0].getCoords().reshape(1,3)
        b_coords = c[1].select(metal_sel)[0].getCoords().reshape(1,3)

        pr.calcTransformation(b_coords, a_coords).apply(c[1])

        outfile = c[0]
        pr.writePDB(workdir + sub_workdir + outfile + '.pdb', c[1])

def get_all_pbd_prody(workdir):
    '''
    find all .pdb file in a folder and load them with prody.
    return a list of pdb_prody.
    '''
    pdbs = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        try:
            pdb_prody = pr.parsePDB(workdir + pdb_path)
            cores.append(pdb_prody)
        except:
            print('not sure')   
    return pdbs


def reduce_dup(pdbs, metal_sel):
    clusters = []
    clustered = set()
    for i in range(len(pdbs)):
        if i in clustered: continue

        cluster = [i]
        clustered.add(i)
        for j in range(i + 1, len(pdbs)):
            if j in clustered: continue
            if len(pdbs[i].select('name CA')) != len(pdbs[j].select('name CA')):
                continue
            print(str(i) + '---' + str(j))
            #TO DO: The calcTransformation here will change the position of pdb. 
            #This will make the output pdb not align well. Current solved by re align.
            pr.calcTransformation(pdbs[j].select('name CA'), pdbs[i].select('name CA')).apply(pdbs[j])
            rmsd = pr.calcRMSD(pdbs[i].select('name CA'), pdbs[j].select('name CA'))

            if rmsd < 0.5:
                cluster.append(j)
                clustered.add(j)

        clusters.append(cluster)
    return clusters

def write_dup_summary(workdir, pdbs, clusters):
    with open(workdir + 'dup_summary.txt', 'w') as f:
        for clu in clusters:
            line = [pdbs[ind].getTitle() for ind in clu]
            f.write('\t'.join(line) + '\n')

def extract_rep_and_writepdb(workdir, pdbs, clusters, sub_workdir = 'reps/'):
    '''
    select represent pdb of each cluster and copy into a subfolder.
    TO DO: It is more efficient to copy the file into the subfolder.
    '''
    if not os.path.exists(workdir + sub_workdir):
        os.mkdir(workdir + sub_workdir)
    sel_pdbs = []

    for clu in clusters:
        sel_pdbs.append((pdbs[clu[0]].getTitle(), pdbs[clu[0]]))

    superimpose_core_and_writepdb(sel_pdbs, sel_pdbs[0], workdir, metal_sel, sub_workdir)

def get_metal_core(file_path, metal_sel):
    try:
        pdb_prody = pr.parsePDB(file_path)
    except:
        print('not sure')
        return []

    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return []
    
    metal_cores = []
    count = 0
    for ni in nis:
        ni_index = ni.getIndex()
        all_near = pdb_prody.select('not water and within 2.6 of index ' + str(ni_index))
        inds = all_near.getResindices()
        count += 1
        metal_cores.append((os.path.basename(file_path) + str(count), pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in inds]))))
            
    return metal_cores

def extract_all_core(workdir, metal_sel):
    '''
    Used for 'all_core_seq' pdbs. 
    all_core_seq pdbs can be produced by extract_all_core_seq()
    '''
    cores = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        core = get_metal_core(workdir + pdb_path, metal_sel)
        cores.extend(core)
    return cores

def writepdb(cores, outdir):
    if not os.path.exists(outdir):
            os.mkdir(outdir)
    for c in cores:
        outfile = c[0]
        pr.writePDB(outdir + outfile + '.pdb', c[1])

def get_aa_core(pdb_prody, aa = 'resname HIS'):
    ni = pdb_prody.select(metal_sel)[0]

    # A pdb can contain more than one NI.
    if not ni:
        return
    
    aa_cores = []
    count = 0

    all_aa = pdb_prody.select(aa)
    if not all_aa:
        return          
    inds = np.unique(all_aa.getResindices())
    for ind in inds:
        count += 1
        sel_pdb_prody = pdb_prody.select('resindex ' + str(ind) + ' '+ str(ni.getResindex()))
        aa_cores.append((pdb_prody.getTitle() + '_' + str(count), sel_pdb_prody))                
    return aa_cores
        
def extract_all_core_aa(pdbs, aa = 'resname HIS'):
    all_aa_cores = []
    for pdb in pdbs:
        aa_cores = get_aa_core(pdb, aa)
        if aa_cores:
            all_aa_cores.extend(aa_cores)
    return all_aa_cores

def superimpose_aa_core(pdbs, outdir, rmsd = 0.5):
    '''
    There are so many ways to superimpose aa and metal.
    This method try to algin the C-CA-N_NI
    '''
    clu = cluster.Cluster()
    clu.rmsd_cutoff = rmsd
    clu.pdbs = []

    for pdb in pdbs:
        c = pdb.select('name C CA N NI').getCoords()
        clu.pdb_coords.append(c)
        clu.pdbs.append(pdb)
    clu.pdb_coords = np.array(clu.pdb_coords, dtype = 'float32')

    clu.make_pairwise_rmsd_mat()  
    if not clu._square:
        clu.make_square()
    if not clu._adj_mat:
        clu.make_adj_mat()
    clu.cluster(min_cluster_size = 2)

    for i in range(len(clu.mems)):
        cluster_out_dir = outdir + str(i) + '/'
        # if not os.path.exists(cluster_out_dir):
        #     os.mkdir(cluster_out_dir)
        print_cluster_pdbs(clu, i, cluster_out_dir, str(rmsd))


def print_cluster_pdbs(clu, rank, outdir='./', tag=''):
    try: os.makedirs(outdir)
    except: pass
    try:
        cent = clu.cents[rank]
        mems = clu.mems[rank]

        # Align backbone of cluster centroid to backbone of centroid of largest cluster.
        R, m_com, t_com = transformation.get_rot_trans(clu.pdb_coords[cent],
                                        clu.pdb_coords[clu.cents[0]])
        cent_coords = np.dot((clu.pdb_coords[cent] - m_com), R) + t_com

        for i, mem in enumerate(mems):
            R, m_com, t_com = transformation.get_rot_trans(clu.pdb_coords[mem], cent_coords)
            pdb = clu.pdbs[mem].copy()
            pdb_coords = pdb.getCoords()
            coords_transformed = np.dot((pdb_coords - m_com), R) + t_com
            pdb.setCoords(coords_transformed)
            is_cent = '_centroid' if mem == cent else ''
            pr.writePDB(outdir + 'cluster_' + str(rank) + '_mem_' + str(i)
                         + is_cent + '_' + pdb.getTitle().split('.')[0] + '.pdb', pdb)

    except IndexError:
        print('Cluster', rank, 'does not exist.')



# work with ProtCID database
'''
import os
import ligand_database as ld

workdir = "/mnt/e/DesignData/ligands/NI/"
metal_sel = 'name NI'

all_pdb_infos = ld.extract_pml_pdb_info(workdir, metal_sel)

#ld.write_pdb_info(workdir + 'protcid_info.txt', all_pdb_infos)

clean_pdb_infos = [info for info in all_pdb_infos if info.selection == 1]

copy_pfam_subfolder(workdir, 'pfam2pdbs_sel/', clean_pdb_infos)


'''
# Compare with rcsb database
'''
pdb_names, pdb_not_in_protcid = ld.compare_rcsb_file(workdir, all_pdb_infos)

workdir_rcsb = '/mnt/e/DesignData/ligands/NI_rcsb/'

pdb_info_not_protcid = ld.extract_rcsb_pdb_info(workdir_rcsb, metal_sel, True, pdb_not_in_protcid)

ld.write_pdb_info(workdir_rcsb + 'pdb_info_not_protcid.txt', pdb_info_not_protcid)

'''

# Extract NI core
''' 
import os
import ligand_database as ld

metal_sel = 'name NI'

workdir = "/mnt/e/DesignData/ligands/NI/pfam2pdbs_sel/"

file_path = workdir + '1a5n4.pfam.pdb'

metal_cores = ld.get_metal_core_seq(file_path, metal_sel)

#ld.extract_all_core(workdir, metal_sel)

ld.extract_all_core_seq(workdir, metal_sel)

'''

# Superimpose (for test purpose)
'''
workdir = "/mnt/e/DesignData/ligands/NI/pfam2pdbs_sel/cores/"


#first = ld.get_metal_core(workdir + '1a5n4.pfam.pdb1_.pdb', metal_sel)[0]

#second = ld.get_metal_core(workdir + '1bxi2.pfam.pdb1_.pdb', metal_sel)[0]

first = pr.parsePDB(workdir + '1a5n4.pfam.pdb1_.pdb')
second = pr.parsePDB(workdir + '1bxi2.pfam.pdb1_.pdb')

first_coords = first.select('name NI')[0].getCoords().reshape(1,3)
second_coords = second.select('name NI')[0].getCoords().reshape(1,3)

pr.calcTransformation(second_coords, first_coords).apply(second)

'''

