import os
import numpy as np
import prody as pr
from scipy.spatial.distance import cdist
from ligand_database import clu_info


def read_cluster_info(file_path):
    clu_infos = []
    if not os.path.exists(file_path):
        return clu_infos
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        for line in lines[1:]:
            ts = line[:-1].split('\t')
            if len(ts)!=7: continue
            clu_infos.append(clu_info(Metal=ts[0], clu_type = ts[1], clu_rmsd= ts[2], total_num=float(ts[3]), clu_rank= int(ts[4]), clu_num= int(ts[5]), score=float(ts[6])))
    return clu_infos
    
def extract_centroid_pdb_in_clu(workdir, file_path = '_summary.txt', score_cut = 0, clu_num_cut = 10):
    pdb_paths = []
    filtered_infos = []
    clu_infos = read_cluster_info(workdir + file_path)
    for info in clu_infos:
        if info.score >= score_cut or info.clu_num >= clu_num_cut:
            filtered_infos.append(info)

    if len(filtered_infos) == 0:
        return pdb_paths

    for info in filtered_infos:
        centroid = [file for file in os.listdir(workdir  + str(info.clu_rank)) if 'centroid' in file][0]
        pdb_paths.append(workdir + str(info.clu_rank) + '/' + centroid)
    return pdb_paths


def _connectivity_filter(arr, inds):
    d = np.diff(arr[inds, :], axis=0)
    tf = (d[:, 0] == 1) & (d[:, 1] == 0) & (d[:, 2] == 0)
    return tf.all()


def connectivity_filter(pdb, window_inds):
    '''
    connectivity_filter copy and modified from qbits.filter
    '''
    N = len(pdb)
    arr = np.zeros((N, 3))
    resnums = pdb.getResnums().astype(np.int16)
    arr[:, 0] = resnums
    chains = pdb.getChids().astype('object')
    arr[:, 1] = np.array(list(map(ord, chains))) 
    segids = pdb.getSegnames().astype('object')
    segids[segids == ''] = 'a' 
    arr[:, 2] = np.array(list(map(ord, segids))) 
    return window_inds[[_connectivity_filter(arr, inds) for inds in window_inds]]


def supperimpose_target_bb(query, target, rmsd_cut = 1):
    '''
    Two possible way:1. Master search; 2. prody calrmsd.
    query: prody pdb
    target: prody pdb
    '''
    query_len = len(query.select('protein and name CA'))
    target_len = len(target.select('protein and name CA'))

    ind = np.arange(query_len)
    window_inds = np.array([ind + i for i in range(target_len- query_len + 1)])

    window_inds = connectivity_filter(target.select('protein and name CA'), window_inds)

    win_extract = []
    for win in window_inds:
        target_sel = target.select('resindex ' + ' '.join([str(w) for w in win]))
        if len(query.select('name N CA C O')) != len(target_sel.select('name N CA C O')):
            continue
        #TO DO: The calcTransformation here will change the position of pdb. 
        #This will make the output pdb not align well. Current solved by re align.
        pr.calcTransformation(query.select('name N CA C O'), target_sel.select('name N CA C O')).apply(query)
        rmsd = pr.calcRMSD(target_sel.select('name N CA C O'), query.select('name N CA C O'))

        if rmsd < rmsd_cut:
            query_curr = query.copy()
            win_extract.append((query_curr, win))

    return win_extract


def write_query_pdbs(outdir, win_extract):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    count = 1
    for wx in win_extract:
        pr.writePDB(outdir + 'ex_' + str(count) + '_' + wx[0].getTitle(), wx[0])
        count += 1

def simple_clash(query, query_2nd):
    '''
    If the two query has CA within 2, then it is a crash.
    '''
    xyzs = []
    for c in query.select('protein and heavy').getCoords():
        xyzs.append(c)
    for c in query_2nd.select('protein and heavy').getCoords():
        xyzs.append(c)

    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    np.fill_diagonal(dists, 5)
    extracts = np.argwhere(dists <= 1.2)

    first_len = len(query.select('protein and heavy'))
    extracts = [(ex[0], ex[1] - first_len) for ex in extracts if ex[0] < first_len and ex[1] > first_len]
    if len(extracts) > 0:
        return True
    return False

def query_target_clash(query, target):
    xyzs = []
    for c in query.select('sc and heavy').getCoords():
        xyzs.append(c)
    for c in target.select('bb').getCoords():
        xyzs.append(c)

    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    np.fill_diagonal(dists, 5)
    extracts = np.argwhere(dists <= 1.2)

    first_len = len(query.select('sc and heavy'))
    extracts = [(ex[0], ex[1] - first_len) for ex in extracts if ex[0] < first_len and ex[1] > first_len]
    if len(extracts) > 0:
        return True
    return False

def geometry_filter():
    '''
    There are only certain geometries for metal contact atoms.
    '''

def metal_distance_extract(target, win_extract, win_extract_2nd, distance_cut = 1):
    xyzs = []
    for win in win_extract:
        xyzs.append(win[0].select('ion').getCoords())
    for win in win_extract_2nd:
        xyzs.append(win[0].select('ion').getCoords())

    xyzs = np.vstack(xyzs)  
    dists = cdist(xyzs, xyzs)

    np.fill_diagonal(dists, 5)
    extracts = np.argwhere(dists <= 1)

    extracts = [(ex[0], ex[1] - len(win_extract)) for ex in extracts if ex[0] < len(win_extract) and ex[1] > len(win_extract)]
    
    extracts_filtered = []
    for i, j in extracts:
        if simple_clash(win_extract[i][0], win_extract_2nd[j][0]): 
            print('two query clash.')
            continue
        if query_target_clash(win_extract[i][0], target) or query_target_clash(win_extract_2nd[j][0], target) :
            print('query target clash.')
            continue
        extracts_filtered.append((i, j))

    return extracts_filtered

def write_cores(outdir, win_extract, win_extract_2nd, extracts):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    count = 1

    for i, j in extracts:
        pr.writePDB(outdir + 'core_' + str(count) + '_' + win_extract[i][0].getTitle(), win_extract[i][0])
        pr.writePDB(outdir + 'core_' + str(count) + '_' + win_extract_2nd[j][0].getTitle(), win_extract_2nd[j][0])
        count += 1
    

def build_core(outdir, target, query, query_2nd, rmsd_cut, rmsd_cut_2nd, distance_cut):
    
    print(query.getTitle())

    print(query_2nd.getTitle())
    
    win_extract = supperimpose_target_bb(query, target, rmsd_cut)

    win_extract_2nd = supperimpose_target_bb(query_2nd, target, rmsd_cut_2nd)

    print(len(win_extract))

    print(len(win_extract_2nd))

    if len(win_extract) == 0 or len(win_extract_2nd) ==0: 
        return len(win_extract), len(win_extract_2nd), None

    extracts = metal_distance_extract(target, win_extract, win_extract_2nd, distance_cut)
    if len(extracts) > 0:
        print(outdir)
        write_cores(outdir, win_extract, win_extract_2nd, extracts)

    return len(win_extract), len(win_extract_2nd), extracts



