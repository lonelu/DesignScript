import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
sys.path.append(r'/mnt/e/GitHub_Design/Qbits')
from qbits import cluster, transformation
import numpy as np
import matplotlib.pyplot as plt

@dataclass
class clu_info:
    Metal: str
    clu_type: str
    clu_rmsd: float
    total_num: int
    clu_rank: int
    clu_num: int
    score: float

    def to_tab_string(self):
        clu_info = self.Metal + '\t' + self.clu_type + '\t' + str(self.clu_rmsd) + '\t' + str(self.total_num) + '\t' + str(self.clu_rank) + '\t' + str(self.clu_num) + '\t'+ str(self.score) 
        return clu_info

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

# download rcsb pdb files

def download_pdb(workdir, filename, resolution = 2.5):
    workdir = "/mnt/e/DesignData/ligands/NI_rcsb/"
    all_pdbs = []
    with open(workdir + filename, 'r') as f: 
        for r in f.readline().split('/t'):
            if r[0] == 'Title' or r[0] == '' or r[4]== '' or (',' in r[4]) or float(r[4]) > 2.5:
                continue
            all_pdbs.append(r[0])
    pr.pathPDBFolder(workdir)
    for p in all_pdbs:
        pr.fetchPDBviaFTP(p, compressed = False) 

    # #Then unzip them in linux with:  
    # cd "/mnt/e/DesignData/ligands/NI_rcsb/"
    # gunzip *.gz

# Basic function. 

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
            pdbs.append(pdb_prody)
        except:
            print('not sure')   
    return pdbs

def writepdb(cores, outdir):
    '''
    cores = list of (pdb name, pdb_prody)
    '''

    if not os.path.exists(outdir):
            os.mkdir(outdir)
    for c in cores:
        outfile = c[0].split('.')[0]
        pr.writePDB(outdir + outfile + '.pdb', c[1])

# Write Summary

def write_clu_info(filename, clu_infos):
    '''
    Write information of cluters.
    @ clu_infos: [clu_info]
    '''
    with open(filename, 'w') as f:
        f.write('Metal\tclu_type\tclu_rmsd\ttotal_num\tclust_rank\tclust_num\tscore\n')
        for r in clu_infos:
            f.write(r.to_tab_string() + '\n')  

def plot_clu_info(clu_infos, outplot):
    fig, ax =plt.subplots(1, 1, figsize=(6,4))
    
    x = list(range(1, len(clu_infos) + 1))
    y = [c.score for c in clu_infos]
    ax.plot(x, y)
    ax.hlines(y=0, xmin = 0, xmax = x[-1], color='r')
    ax.legend()
    #ax.set_xticks(x)
    
    ax.set_xlabel("Rank", fontsize = 12)
    ax.set_ylabel("vdM score", fontsize = 12)
    plt.tight_layout()
    plt.savefig(outplot)
    plt.close()
# Prepare rcsb database. // extract seq within +-3 aa for each contact aa. 

def connectivity_filter(pdb_prody, ind, ext_ind):
    res1 = pdb_prody.select('protein and resindex ' + str(ind))
    res2 = pdb_prody.select('protein and resindex ' + str(ext_ind))
    if not res2:
        return False
    if res1[0].getResnum() - res2[0].getResnum() == ind - ext_ind and res1[0].getChid() == res2[0].getChid() and res1[0].getSegname() == res2[0].getSegname():
        return True
    return False

def extend_res_indices(inds_near_res, pdb_prody, extend = 3):
    extend_inds = []
    inds = set()
    for ind in inds_near_res:
        for i in range(-extend, extend + 1):
            the_ind = ind + i
            if the_ind not in inds and the_ind>= 0 and connectivity_filter(pdb_prody, ind, the_ind):
                extend_inds.append(the_ind)
                inds.add(the_ind)
    return extend_inds

def get_metal_core_seq(pdb_prody, metal_sel, extend = 3):
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return
    
    metal_cores = []
    count = 0
    for ni in nis:
        ni_index = ni.getIndex()
        #all_near = pdb_prody.select('nitrogen or oxygen').select('not water and within 2.6 of index ' + str(ni_index))
        all_near = pdb_prody.select('nitrogen or oxygen').select('protein and within 2.6 of index ' + str(ni_index))
        if not all_near or len(all_near) < 3:
            continue          
        inds = all_near.getResindices()
        all_near_res = pdb_prody.select('protein and resindex ' + ' '.join([str(ind) for ind in inds]))
        if not all_near_res or len(np.unique(all_near_res.getResindices())) < 2:
            continue     
        inds_near_res =  all_near_res.getResindices()
        ext_inds = extend_res_indices(inds_near_res, pdb_prody, extend)
        count += 1
        sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
        metal_cores.append((pdb_prody.getTitle() + '_' + str(count), sel_pdb_prody))        
    return metal_cores

def extract_all_core_seq(pdbs, metal_sel, extend = 3):
    cores = []
    for pdb in pdbs:
        core = get_metal_core_seq(pdb, metal_sel, extend)
        if core:
            cores.extend(core)
    return cores

def superimpose_core_and_writepdb(cores, first, metal_sel, outdir):
    '''
    Just supperimpose on the metal. And write pdb.
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for c in cores:
        a_coords = first[1].select(metal_sel)[0].getCoords().reshape(1,3)
        b_coords = c[1].select(metal_sel)[0].getCoords().reshape(1,3)

        pr.calcTransformation(b_coords, a_coords).apply(c[1])

        outfile = c[0]
        pr.writePDB(outdir + outfile + '.pdb', c[1])

# Prepare rcsb database. // Reduce duplication.

def reduce_dup(pdbs, metal_sel):
    '''
    Cluster pdbs based on the structure similarity.
    '''
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
            pr.calcTransformation(pdbs[j].select('name N CA C O'), pdbs[i].select('name N CA C O')).apply(pdbs[j])
            rmsd = pr.calcRMSD(pdbs[i].select('name N CA C O'), pdbs[j].select('name N CA C O'))

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

def extract_rep_and_writepdb(pdbs, clusters, metal_sel, outdir):
    '''
    select represent pdb of each cluster and copy into a subfolder.
    TO DO: It is more efficient to copy the file into the subfolder.
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    sel_pdbs = []

    for clu in clusters:
        sel_pdbs.append((pdbs[clu[0]].getTitle(), pdbs[clu[0]]))

    superimpose_core_and_writepdb(sel_pdbs, sel_pdbs[0], metal_sel, outdir)

# Extract Cores from the prepared database. //For manucheck. 

def get_metal_core(pdb_prody, metal_sel):
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
        metal_cores.append((pdb_prody.getTitle(), pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in inds]))))
            
    return metal_cores

def extract_all_core(pdbs, metal_sel):
    cores = []
    for pdb in pdbs:
        core = get_metal_core(pdb, metal_sel)
        if core:
            cores.extend(core)
    return cores

# Extract and Cluster Cores from the prepared database.

def get_aa_core(pdb_prody, metal_sel, aa = 'resname HIS', consider_phipsi = False, extention = 0):
    '''
    extract amino acid core + metal.
    if consider_phipsi == True. Extract 
    '''
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        print('No NI?' + pdb_prody.getTitle())
        return
    
    ni = nis[0]
    aa_cores = []
    count = 0

    all_aa = pdb_prody.select(aa + ' and within 2.6 of index ' + str(ni.getIndex()))
    if not all_aa:
        return          
    inds = np.unique(all_aa.getResindices())
    for ind in inds:
        count += 1
        if consider_phipsi:
            ext_inds = extend_res_indices([ind], pdb_prody, extend =1)
            if len(ext_inds) != 3:
                continue
            print(pdb_prody.getTitle() + '+' + '-'.join([str(x) for x in ext_inds]))
            atom_inds = []
            atom_inds.extend(pdb_prody.select('resindex ' + str(ind-1)).select('name C O').getIndices())
            atom_inds.extend(pdb_prody.select('resindex ' + str(ind)).getIndices())
            atom_inds.extend(pdb_prody.select('resindex ' + str(ind+1)).select('name N CA').getIndices())
            sel_pdb_prody = pdb_prody.select('index ' + ' '.join([str(x) for x in atom_inds]) + ' '+ str(ni.getIndex()))
        if extention != 0:
            ext_inds = extend_res_indices([ind], pdb_prody, extend =extention)
            if len(ext_inds) != 2*extention + 1:
                continue
            print(pdb_prody.getTitle() + '+' + '-'.join([str(x) for x in ext_inds]))
            sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
        else:
            sel_pdb_prody = pdb_prody.select('resindex ' + str(ind) + ' '+ str(ni.getResindex()))
        aa_cores.append((pdb_prody.getTitle() + '_' + str(count), sel_pdb_prody))                
    return aa_cores
        
def extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', consider_phipsi = False, extention = 0):
    all_aa_cores = []
    for pdb in pdbs:
        aa_cores = get_aa_core(pdb, metal_sel, aa, consider_phipsi, extention)

        if aa_cores:
            all_aa_cores.extend(aa_cores)
    return all_aa_cores

def superimpose_aa_core(pdbs, rmsd = 0.5, len_sel = 5, align_sel = 'name C CA N O NI'):
    '''
    There are so many ways to superimpose aa and metal.
    This method try to algin the C-CA-N_NI
    '''
    clu = cluster.Cluster()
    clu.rmsd_cutoff = rmsd
    clu.pdbs = []

    for pdb in pdbs:
        c = pdb.select(align_sel).getCoords()       
        if len(c)!= len_sel: 
            continue
        clu.pdb_coords.append(c)
        clu.pdbs.append(pdb)
    clu.pdb_coords = np.array(clu.pdb_coords, dtype = 'float32')

    clu.make_pairwise_rmsd_mat()  
    if not clu._square:
        clu.make_square()
    if not clu._adj_mat:
        clu.make_adj_mat()
    clu.cluster(min_cluster_size = 2)

    return clu

def get_clu_info_write(outfile, pdbs, clu, rmsd = 0.5, metal_sel = 'NI', align_sel = 'name C CA N O NI'):
    clu_infos = []
    n_avg = sum([len(m) for m in clu.mems])/len(clu.mems)
    for i in range(len(clu.mems)):
        c = clu_info(Metal=metal_sel, clu_type = align_sel, 
            clu_rmsd=rmsd, total_num = len(pdbs), clu_rank = i, 
            clu_num=len(clu.mems[i]), score = np.log(len(clu.mems[i])/n_avg) )
        clu_infos.append(c)
    write_clu_info(outfile, clu_infos)
    return clu_infos

def print_cluster_pdbs(clu, outdir, rmsd = 0.5):
    for i in range(len(clu.mems)):
        cluster_out_dir = outdir + str(i) + '/'
        # if not os.path.exists(cluster_out_dir):
        #     os.mkdir(cluster_out_dir)
        _print_cluster_rank_pdbs(clu, i, cluster_out_dir, str(rmsd))

def _print_cluster_rank_pdbs(clu, rank, outdir='./', tag=''):
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

       