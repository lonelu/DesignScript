import os
import prody as pr
from scipy.spatial.distance import cdist
from dataclasses import dataclass
import shutil
import sys
import numpy as np
import matplotlib.pyplot as plt
import cluster

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
                all_lines.extend(f.readlines())
    with open(workdir + 'all_rcsb.txt', 'w') as f:
        f.write('\t'.join(all_lines[0].split(',')))
        for r in all_lines:
            if 'Entry ID' not in r and r.split(',')[0]!= '':
                f.write('\t'.join(r.split(',')))

# download rcsb pdb files

def download_pdb(workdir, filename, resolution = 2.5):
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    all_pdbs = []
    with open(workdir + filename, 'r') as f: 
        for line in f.readlines():
            #print(line)
            r = line.split('\t')
            #print(r)
            if r[0] == '"Entry ID"': continue
            if r[0] == '' or r[4]== '' or (',' in r[4]) or float(r[4].split('"')[1]) > 2.5:
                continue
            all_pdbs.append(r[0].split('"')[1])

    exist_pdb = set()
    for file in os.listdir(workdir):
        if file.endswith(".pdb.gz"):
            exist_pdb.add(file.split('.')[0].upper())

    pr.pathPDBFolder(workdir)
    for p in all_pdbs:
        if p in exist_pdb: continue
        pr.fetchPDBviaFTP(p, compressed = False) 

    # #Then unzip them in linux with:  
    # cd /mnt/e/DesignData/ligands/NI_rcsb/
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
    fig, (ax, ax1) =plt.subplots(2, 1, figsize=(12,8))
    
    x = list(range(1, len(clu_infos) + 1))
    y = [c.score for c in clu_infos]
    ax.plot(x, y)
    ax.hlines(y=0, xmin = 0, xmax = x[-1], color='r')
    ax.legend()
    #ax.set_xticks(x)
    ax.set_xlabel("Rank", fontsize = 12)
    ax.set_ylabel("vdM score", fontsize = 12)

    counts = [c.clu_num for c in clu_infos]
    ax1.bar(x, counts)

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

def extend_res_indices(inds_near_res, pdb_prody, extend = 4):
    extend_inds = []
    inds = set()
    for ind in inds_near_res:
        for i in range(-extend, extend + 1):
            the_ind = ind + i
            if the_ind not in inds and the_ind>= 0 and connectivity_filter(pdb_prody, ind, the_ind):
                extend_inds.append(the_ind)
                inds.add(the_ind)
    return extend_inds

def get_metal_core_seq(pdb_prody, metal_sel, extend = 4):
    metal = metal_sel.split(' ')[-1]
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        return
    
    metal_cores = []
    count = 0
    for ni in nis:
        ni_index = ni.getIndex()
        #all_near = pdb_prody.select('nitrogen or oxygen or sulfur').select('not water and within 2.83 of index ' + str(ni_index))
        all_near = pdb_prody.select('protein and within 2.83 of index ' + str(ni_index))
        if not all_near or not all_near.select('nitrogen or oxygen or sulfur') or len(all_near.select('nitrogen or oxygen or sulfur')) < 3:
            continue          
        inds = all_near.select('nitrogen or oxygen or sulfur').getResindices()
        # all_near_res = pdb_prody.select('protein and resindex ' + ' '.join([str(ind) for ind in inds]))
        # if not all_near_res or len(np.unique(all_near_res.getResindices())) < 2:
        #     continue     
        # inds_near_res =  all_near_res.getResindices()
        # ext_inds = extend_res_indices(inds_near_res, pdb_prody, extend)
        ext_inds = extend_res_indices(inds, pdb_prody, extend)
        count += 1
        sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
        metal_cores.append((pdb_prody.getTitle() + '_' + metal + '_'+ str(count), sel_pdb_prody))        
    return metal_cores

def extract_all_core_seq_from_path(workdir, metal_sel, extend = 3):
    cores = []
    for pdb_path in os.listdir(workdir):
        if not pdb_path.endswith(".pdb"):
            continue
        try:
            pdb_prody = pr.parsePDB(workdir + pdb_path)
            core = get_metal_core_seq(pdb_prody, metal_sel, extend)
            if core:
                cores.extend(core)
            pdb_prody = None
        except:
            print('not sure')   
    return cores

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
            if len(pdbs[i].select('name N CA C O')) != len(pdbs[j].select('name N CA C O')):
                continue
            #print(str(i) + '---' + str(j))
            sequence_equal = False
            if pdbs[i].select('name CA').getSequence() == pdbs[i].select('name CA').getSequence():
                sequence_equal = True
            #TO DO: The calcTransformation here will change the position of pdb. 
            #This will make the output pdb not align well. Current solved by re align.
            pr.calcTransformation(pdbs[j].select('name N CA C O'), pdbs[i].select('name N CA C O')).apply(pdbs[j])
            rmsd = pr.calcRMSD(pdbs[i].select('name N CA C O'), pdbs[j].select('name N CA C O'))

            if (sequence_equal and rmsd < 0.8) or rmsd < 0.25:
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
        all_near = pdb_prody.select('not water and within 2.83 of index ' + str(ni_index))
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
    aa_name = aa.split(' ')[-1]
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        print('No NI?' + pdb_prody.getTitle())
        return
    
    ni = nis[0]
    aa_cores = []
    count = 0

    all_aa = pdb_prody.select(aa + ' and within 2.83 of index ' + str(ni.getIndex()))
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
        elif extention != 0:
            ext_inds = extend_res_indices([ind], pdb_prody, extend =extention)
            if len(ext_inds) != 2*extention + 1:
                continue
            print(pdb_prody.getTitle() + '+' + '-'.join([str(x) for x in ext_inds]))
            sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
        else:
            sel_pdb_prody = pdb_prody.select('resindex ' + str(ind) + ' '+ str(ni.getResindex()))
        aa_cores.append((pdb_prody.getTitle() + '_' + aa_name + '_' + str(count), sel_pdb_prody))                
    return aa_cores
        
def get_2aa_core(pdb_prody, metal_sel, aa = 'resname HIS', extention = 3, extention_out = 1):
    '''
    extract 2 amino acid core + metal.

    '''
    aa_name = aa.split(' ')[-1]
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        print('No NI?' + pdb_prody.getTitle())
        return
    
    ni = nis[0]
    aa_cores = []
    count = 0

    all_aa = pdb_prody.select(aa + ' and within 2.83 of index ' + str(ni.getIndex()))
    if not all_aa:
        return          
    inds = np.unique(all_aa.getResindices())
    exts = []
    for ind in inds:
        ext_inds = extend_res_indices([ind], pdb_prody, extend =extention)
        exts.append(ext_inds)
    
    pairs = []
    pairs_ind = []
    for i in range(len(inds) - 1):
        for j in range(i+1, len(inds)):
            overlap = list(set(exts[i]) & set(exts[j]))
            if len(overlap) ==0: continue
            pairs.append((i, j))
            if inds[i] < inds[j]:
                pairs_ind.append((inds[i], inds[j]))
            else:
                pairs_ind.append((inds[j], inds[i]))
    
    for v in range(len(pairs)):
        i = pairs[v][0]
        j = pairs[v][1]
        count += 1
        ext_inds = overlap = list(set(exts[i]) | set(exts[j]))
        ext_inds = [x for x in ext_inds if x >= pairs_ind[v][0]-extention_out and x <= pairs_ind[v][1]+extention_out]
        sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + ' '+ str(ni.getResindex()))
        aa_cores.append((pdb_prody.getTitle() + '_' + aa_name + '_' + str(count), sel_pdb_prody))                
    return aa_cores

def extract_all_core_aa(pdbs, metal_sel, aa = 'resname HIS', consider_phipsi = False, extention = 0, extract2aa = False):
    all_aa_cores = []
    for pdb in pdbs:
        if extract2aa:
            aa_cores = get_2aa_core(pdb, metal_sel, aa, extention)
        else:
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
    if len(clu.pdb_coords) <= 0:
        return 
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

def check_metal_diff(clu, metals = ['CA', 'CO', 'CU', 'FE', 'MG', 'MN', 'NI', 'ZN']):
    metal_diffs = []
    for rank in range(len(clu.mems)):
        if not clu.mems[rank].any(): continue

        mems = clu.mems[rank]  
        metal_diff = [0]*len(metals)
        for i, mem in enumerate(mems):
            pdb = clu.pdbs[mem]

            if pdb.select('ion and name CA'):
                metal_diff[0]+=1
            elif pdb.select('name CO'):
                metal_diff[1]+=1
            elif pdb.select('name CU'):
                metal_diff[2]+=1
            elif pdb.select('name FE'):
                metal_diff[3]+=1
            elif pdb.select('name MG'):
                metal_diff[4]+=1                
            elif pdb.select('name MN'):
                metal_diff[5]+=1
            elif pdb.select('name NI'):
                metal_diff[6]+=1
            elif pdb.select('name ZN'):
                metal_diff[7]+=1

        metal_diffs.append(metal_diff)
    return metal_diffs    

def write_metal_diff(workdir, metal_diffs, metals = ['CA', 'CO', 'CU', 'FE', 'MG', 'MN', 'NI', 'ZN']):
    with open(workdir + 'metal_diff.txt', 'w') as f:
        f.write('\t'.join(metals) + '\n')
        for md in metal_diffs:
            f.write('\t'.join([str(x) for x in md]) + '\n')

# Extract 3d angle

def get_atom_core(pdb_prody, metal_sel, tag = '_atom'):
    '''
    extract all nearby atoms
    '''
    nis = pdb_prody.select(metal_sel)

    # A pdb can contain more than one NI.
    if not nis:
        print('No Metal?' + pdb_prody.getTitle())
        return
    
    ni = nis[0]
    atom_cores = []
    count = 0

    all_atom = pdb_prody.select('protein and within 2.83 of index ' + str(ni.getIndex())).select('nitrogen or oxygen or sulfur')
    if not all_atom:
        return     
    atom_inds = np.unique(all_atom.getIndices())
    sel_pdb_prody = pdb_prody.select('index ' + ' '.join([str(x) for x in atom_inds]) + ' '+ str(ni.getIndex()))
    atom_cores.append(( pdb_prody.getTitle() + tag, sel_pdb_prody))   

    return atom_cores

def extract_all_atom_core(pdbs, metal_sel, tag  ='_atom'):
    all_atom_cores = []
    for pdb in pdbs:
        atom_cores = get_atom_core(pdb, metal_sel, tag)
        if atom_cores:
            all_atom_cores.extend(atom_cores)
    return all_atom_cores


# run cluster

def run_cluster(_pdbs, workdir, outdir, rmsd, metal_sel, len_sel, align_sel):
    
    clu = superimpose_aa_core(_pdbs, rmsd = rmsd, len_sel = len_sel, align_sel = align_sel)
    
    if not clu or len(clu.mems) == 0: return
    
    print_cluster_pdbs(clu, workdir + outdir)

    metal_diffs = check_metal_diff(clu)

    write_metal_diff(workdir + outdir, metal_diffs)

    clu_infos = get_clu_info_write(workdir + outdir + '_summary.txt', _pdbs, clu, rmsd = rmsd, metal_sel = metal_sel, align_sel = align_sel)

    plot_clu_info(clu_infos, workdir + outdir + '_score.png')