import prody as pr
import numpy as np
import pickle
import os
import sys
sys.path.append(r'/mnt/e/GitHub_Design/DesignScript/topology/')
import topology

'''
python /mnt/e/GitHub_Design/DesignScript/topology/topology_usage.py
'''

def rec_dd():
    """returns a recursive dictionary"""
    return topology.defaultdict(rec_dd)
    
def run_topology(workdir, outdir, pdb_file, constrained_rns, constrained_rns_vals, selection, resfile_pre= None):
    os.makedirs(outdir, exist_ok = True)

    pdb_path = workdir + pdb_file
    pdb = pr.parsePDB(pdb_path)
    #>>>**list of constrained residues paired with list of res type (-1, 0, or 1)
    top = topology.Topology(**dict(constrained_rns=constrained_rns, constrained_rns_vals=constrained_rns_vals))
    #>>> list of residues per helix ***important that helix 2 and helix 4 are same length**

    top.load_pdb(pdb, selection=selection)
    top.load_pdb_ala(pdb, selection=selection)
    top.set_topologies(outdir=workdir)
    top.set_surface_res()
    top.set_contacts()
    top.run_mc(num_iterations=400000)
    seqs = [top.seqs[i] for i in range(len(top.seqs)) if np.sum(list(top.seqs[i].values())) <= -2]
    en_gaps = [top.en_gaps[i] for i in range(len(top.seqs)) if np.sum(list(top.seqs[i].values())) <= -2]
    en_fs = [top.en_fs[i] for i in range(len(top.seqs)) if np.sum(list(top.seqs[i].values())) <= -2]
    seq_rep_lens = [top.seq_rep_lens[i] for i in range(len(top.seqs)) if np.sum(list(top.seqs[i].values())) <= -2]
    top.seq_rep_lens = seq_rep_lens
    top.en_gaps = en_gaps
    top.en_fs = en_fs
    top.seqs = seqs
    top.find_pareto_front()
    top.find_nearest_utopian_pt(weight_en_f=0.75, weight_seq_rep_len=0.5)
    top.map_seq_resnums_to_pdb(pdb)
    top.set_charge_groups(top.seqs[top.nearest_utopian_pt])
    top.print_charge_groups()
    top.save_sequence(top.seqs[top.nearest_utopian_pt], outdir=outdir)
    pr.execDSSP(pdb_path, outputdir = outdir)


    dssp = topology.parse_dssp(outdir + pdb_file.split('.')[0] + '.dssp')

    with open('/mnt/e/GitHub_Design/DesignScript/topology/ss_burial_propensity_label_dict_20180902.pkl', 'rb') as infile:
        ss_bur_prop_label_dict = pickle.load(infile)
    polar_dict = {}
    polar_dict[-1] = set('ED')
    polar_dict[1] = set('KR')
    polar_dict[0] = polar_dict[-1] | polar_dict[1]

    with open(outdir + 'surface_sequence.txt','r') as infile:
        infile.readline()
        surface_seq = {}
        for line in infile:
            line = line.split()
            surface_seq[int(line[0])] = int(line[1])

    resfile = {}
    if resfile_pre:
        with open('resfile.txt', 'r') as infile:
            infile.readline()
            for line in infile:
                try:
                    line = line.strip().split()
                    resfile[(int(line[0]), line[1])] = line[2:]
                except:
                    pass
    pdb_ala = pr.parsePDB(pdb_path)
    resind_exp, resind_int, resind_bur = topology.partition_res_by_burial(pdb_ala, alpha=9)
    resnum_exp = set(pdb_ala.select('resindex ' + ' '.join([str(r) for r in resind_exp])).getResnums())
    resnum_int = set(pdb_ala.select('resindex ' + ' '.join([str(r) for r in resind_int])).getResnums())
    resnum_bur = set(pdb_ala.select('resindex ' + ' '.join([str(r) for r in resind_bur])).getResnums())
    dontallow_global = set('PCHM')
    add_to_all = set('AV') #set('AV')
    res_burial_map = {}
    for res in resnum_exp:
        res_burial_map[res] = 'e'
    for res in resnum_int:
        res_burial_map[res] = 'i'
    for res in resnum_bur:
        res_burial_map[res] = 'b'
    with open(outdir + 'resfile_topology.txt', 'a') as outfile:
        for (res, ch), ss in dssp.items():
            if (res, ch) in resfile.keys():
                continue
            bur = res_burial_map[res]
            aa_set = {aa for aa, prop in ss_bur_prop_label_dict[ss][bur].items() if np.round(prop, 1) >= 0.9} - dontallow_global
            aa_set |= add_to_all
            if res in surface_seq.keys():
                q = surface_seq[res]
                aa_set = aa_set - polar_dict[-1*q]
            line = str(res) + ' ' + ch + ' PIKAA ' + ''.join(aa_set) + ' \n'
            outfile.write(line)
    return

#>>>
workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/search_unknow_result/output__helix6a_10-14-139_HDH.pdb_20220624-222024/loop_ss_20220625-103345/merge/'
outdir = workdir + 'output/'
pdb_file = 'helix6a_10-14-139_tts_rosetta_prep_all_ala_res.pdb'
constrained_rns=[10, 11, 14, 58, 76, 79, 116]
constrained_rns_vals=[0, 0, 0, 0, 0, 0, 0]
selection='resnum 6to26 36to56 71to91 101to121'

#>>>
workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/design_0/'
outdir = workdir + 'output/'
pdb_file = 'helix6_loop1_cut.pdb'
constrained_rns=[]
constrained_rns_vals=[]
selection='resnum 4to34 41to71 79to109 115to145'

#>>>
workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/design_1/'
outdir = workdir + 'output/'
pdb_file = 'helix6_loop2_cut.pdb'
constrained_rns=[]
constrained_rns_vals=[]
selection='resnum 3to26 33to56 64to87 95to118'

#>>>
workdir = '/mnt/e/DesignData/_temp/yang/20220726/'
outdir = workdir + 'output/'
pdb_file = '6w70a_ala3.pdb'
constrained_rns=[]
constrained_rns_vals=[]
selection='resnum 3to26 33to56 64to87 95to118'

#>>>
workdir = '/mnt/e/DesignData/Metalloenzyme/SAHA_Vorinostat/run_design_cgs3/SAHA_Rosetta20220728/Rosetta/'
outdir = workdir + 'output/'
pdb_file = 'bb_prep_A.pdb'
constrained_rns=[]
constrained_rns_vals=[]
selection='resnum 4to32 43to71 80to108 117to145'

#>>> C3 6 helix bundle is not working
workdir = '/mnt/e/DesignData/tm/Kehan/20220804/_design/'
outdir = workdir + 'topo_output/'
pdb_file = 'bb_prep_loop_A.pdb'
constrained_rns=[]
constrained_rns_vals=[]
selection=None

#>>> fluo
workdir = '/mnt/e/DesignData/bpp_fluo_comb/fluo/output1_09_f63440_nick_ala/sel/Rosetta/'
outdir = workdir + 'topo_output1/'
pdb_file = 'bb_prep_ALA_topo.pdb'
constrained_rns=[]
constrained_rns_vals=[]
selection='resnum 5to26 33to54 59to80 87to108'

#>>> 44b
workdir = '/mnt/e/DesignData/bpp_fluo_comb/44b/44b_kp/'
outdir = workdir + 'topo_output/'
pdb_file = 'bb_prep_ALA_topo_A.pdb'
constrained_rns=[]
constrained_rns_vals=[]
selection='resnum 1to19 26to44 52to70 82to100'

#>>> short599
workdir = workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC_short/loop5/Rosetta/'
outdir = workdir + 'topo_output/'
pdb_file = 'bb_prep_ALA_topo.pdb'
constrained_rns=[]
constrained_rns_vals=[]
selection='resnum 5to22 27to44 52to69 73to90'

run_topology(workdir, outdir, pdb_file, constrained_rns, constrained_rns_vals, selection, resfile_pre= None)