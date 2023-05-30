import random
import prody as pr
import numpy as np
from collections import defaultdict, deque
from itertools import groupby
#from .convex_hull import partition_res_by_burial
from scipy.spatial.distance import cdist
import traceback

#from .convex_hull
def partition_res_by_burial(pdb_ala, alpha=9):
    """Returns residue indices of exposed, intermediate, and buried residues
    based on CA hull and CB hull."""
    ahull_ca = AlphaHull(alpha=alpha)
    ahull_ca.coords = pdb_ala.select('name CA').getCoords()
    ahull_ca.calc_hull()
    ahull_cb = AlphaHull(alpha=alpha)
    ahull_cb.set_coords(pdb_ala)
    ahull_cb.calc_hull()    
    ahull_cb.set_resindices(pdb_ala)
    cb_in_ca_hull = ahull_ca.pnts_in_hull(ahull_cb.coords)
    resindices_cb_in_ca_hull = set(ahull_cb.resindices[cb_in_ca_hull])
    resindices_cb_hull = set(ahull_cb.resindices[np.unique(ahull_cb.hull)])
    resindices_not_cb_hull = set(ahull_cb.resindices) - resindices_cb_hull
    resindices_exposed = resindices_cb_hull - resindices_cb_in_ca_hull
    resindices_intermediate = resindices_cb_in_ca_hull - resindices_not_cb_hull
    resindices_buried = resindices_cb_in_ca_hull - resindices_intermediate
    res_ = resindices_not_cb_hull - resindices_buried
    resindices_intermediate |= res_
    return resindices_exposed, resindices_intermediate, resindices_buried

# Energy table for Monte Carlo. -1 is neg, 0 neutral, +1 positive
En = defaultdict(dict)
En[-1][-1] = 3
En[1][1] = 2
En[1][-1] = -1
En[-1][1] = -1
En[1][0] = -0.1
En[0][1] = -0.1
En[0][0] = 0
En[0][-1] = -0.1
En[-1][0] = -0.1


class Topology:
    """Calculates surface charge distribution for helical bundles to
    stabilize forward topology over reverse topology.

    Example of typical usage:

    top = combs.apps.topology.Topology()
    top.load_pdb(pdb, selection='resnum 1to31 39to68 79to106 111to140')
    top.set_topologies(outdir='surf_test/')
    top.set_surface_res()
    top.set_contacts()
    top.run_mc()
    top.find_pareto_front()
    top.map_seq_resnums_to_pdb(top.pdb_f)
    top.set_charge_groups(top.seqs[162063])
    top.print_charge_groups()

    """
    def __init__(self, **kwargs):
        self.pdb = None
        self.pdb_sel = None
        self.pdb_ala = None
        self.pdb_ala_sel = None
        self.pdb_f = None
        self.pdb_r = None
        self.surf_rns_f = None
        self.surf_rns_r = None
        self.surf_rns = None
        self.surf_sel_f = None
        self.surf_sel_r = None
        self.nbrs_f = None
        self.nbrs_r = None
        self.contacts_f = None
        self.contacts_r = None
        self.seq = dict()
        self.seqs = list()
        self.en_gaps = list()
        self.en_fs = list()
        self.seq_rep_lens = list()
        self.resnum_conv = dict()
        self.resnum_conv_reverse = dict()
        self.negs = list()
        self.neuts = list()
        self.poss = list()
        self._surf_rns = None
        self.pareto_front = None
        self.nearest_utopian_pt = None
        self.constrained_rns = kwargs.get('constrained_rns', list())
        self.constrained_rns_vals = kwargs.get('constrained_rns_vals', list())
        self.constrained_surf_rns_f = list()
        self.constrained_surf_rns_vals_f = dict()

    def load_pdb(self, pdb, selection=None):
        """Loads a prody pdb object (pdb). Takes a selection string (selection).
        Selection string is needed if the pdb is not disconnected helices. Note
        that the selection needs to contain swapped helices of equal lengths, e.g.
        helices 1 and 3 of a 4-helix bundle (with helices 0,1,2,3) must be the
        same length"""
        pdb = pdb.select('protein').copy()
        self.pdb = pdb
        if selection is not None:
            self.pdb_sel = pdb.select(selection)
        else:
            self.pdb_sel = pdb

    def load_pdb_ala(self, pdb, selection=None):
        """Loads a prody pdb object (pdb) that is all alanine residues. Takes a
        selection string (selection).
        Selection string is needed if the pdb is not disconnected helices. Note
        that the selection needs to contain swapped helices of equal lengths, e.g.
        helices 1 and 3 of a 4-helix bundle (with helices 0,1,2,3) must be the
        same length"""
        pdb = pdb.select('protein').copy()
        if (len(set(pdb.getResnames())) != 1) or (set(pdb.getResnames()).pop() != 'ALA'):
            raise "*pdb* must be all alanine residues."
        self.pdb_ala = pdb
        if selection is not None:
            self.pdb_ala_sel = pdb.select(selection)
        else:
            self.pdb_ala_sel = pdb

    def set_topologies(self, outdir=None, tag=''):
        """Takes a prody pdb object or prody selection.
        Outputs two prody pdb objects: pdb with topology 1 and
        pdb with topology 2. This function finds the best cyclic permutation
        of the helices so that helix 0 and n (if n exists, e.g. n = 2
        in a 4 helix bundle, n = 3 in a 6 helix bundle) can have arbitrary length, but
        other helices must be the same length (so that they can be swapped
        in the structure).  Ideally, the pdb should have CB atoms for
        subsequent alpha hull calculations."""
        self.pdb_f = number_helices(self.pdb_sel, reverse=False)
        self.pdb_r = number_helices(self.pdb_sel, reverse=True)
        if outdir is not None:
            if outdir[-1] != '/':
                outdir += '/'
            pr.writePDB(outdir + 'pdb_f' + tag + '.pdb', self.pdb_f)
            pr.writePDB(outdir + 'pdb_r' + tag + '.pdb', self.pdb_r)

    @staticmethod
    def _map_resnums_to_pdb(resnums_of_pdb1, pdb1, pdb2):
        """Maps resnums of pdb1 (list of integers, resnums_of_pdb1) to corresponding
        resnums of pdb2. Returns forward and reverse resnum conversion dictionaries."""
        resnum_conv = dict()
        resnum_conv_reverse = dict()
        for rn in resnums_of_pdb1:
            try:
                sel = pdb1.select('name CA and resnum ' + str(rn))
                pdb_rn = pdb2.select('within 0.05 of sel', sel=sel).getResnums()[0]
                resnum_conv[rn] = pdb_rn
                resnum_conv_reverse[pdb_rn] = rn
            except:
                traceback.print_exc()
        return resnum_conv, resnum_conv_reverse

    def set_surface_res(self, alpha=9, selection_only=False):
        """Calculates surface residues of forward and reverse topologies
        based on alpha hull calculation."""
        if selection_only:
            exp, inter, bur = partition_res_by_burial(self.pdb_ala_sel, alpha=alpha)
            rns = self.pdb_ala_sel.select('name CB or (resname GLY and name CA)').getResnums()
        else:
            exp, inter, bur = partition_res_by_burial(self.pdb_ala, alpha=alpha)
            rns = self.pdb_ala.select('name CB or (resname GLY and name CA)').getResnums()

        self._surf_rns = [rns[rn] for rn in exp]
        resnum_conv_f, resnum_conv_reverse_f = self._map_resnums_to_pdb(self._surf_rns,
                                                                        self.pdb, self.pdb_f)
        self.surf_rns_f = set(list(resnum_conv_f.values()))
        resnum_conv_r, resnum_conv_reverse_r = self._map_resnums_to_pdb(self._surf_rns,
                                                                        self.pdb, self.pdb_r)
        self.surf_rns_r = set(list(resnum_conv_r.values()))
        self.surf_rns = list(self.surf_rns_f) # list(self.surf_rns_f | self.surf_rns_r)
        self.surf_sel_f = self.pdb_f.select('name CA and resnum ' +
                                            ' '.join(str(rn) for rn in self.surf_rns))
        self.surf_sel_r = self.pdb_r.select('name CA and resnum ' +
                                            ' '.join(str(rn) for rn in self.surf_rns))

        if self.constrained_rns:
            constrained_surf_rns = [c_rn for c_rn in self.constrained_rns if c_rn in self._surf_rns]
            constrained_surf_rns_vals = [val for c_rn, val in zip(self.constrained_rns, self.constrained_rns_vals)
                                              if c_rn in self._surf_rns]
            if constrained_surf_rns:
                constrained_surf_rns_f, resnum_cons_conv_reverse_f = self._map_resnums_to_pdb(constrained_surf_rns,
                                                                                              self.pdb, self.pdb_f)
                self.constrained_surf_rns_f = list(constrained_surf_rns_f.values())
                for key, val in zip(self.constrained_surf_rns_f, constrained_surf_rns_vals):
                    self.constrained_surf_rns_vals_f[key] = val


    def set_contacts(self, calpha_distance=10.5):
        """Sets contacts between surface residues. A contact is by default defined
        as a pair of C alpha atoms with a distance less than 10 angstroms."""
        self.nbrs_f = pr.findNeighbors(self.surf_sel_f, calpha_distance)
        self.nbrs_r = pr.findNeighbors(self.surf_sel_r, calpha_distance)
        self.contacts_f = list()
        for nbr in self.nbrs_f:
            resnum_0 = nbr[0].getResnum()
            resnum_1 = nbr[1].getResnum()
            filter = True
            if np.abs(resnum_0 - resnum_1) > 6:
                resind_0 = nbr[0].getResindex()
                resind_1 = nbr[1].getResindex()
                ca_0 = nbr[0]
                cb_0 = self.pdb_f.select('name CB and resindex ' + str(resind_0))
                ca_1 = nbr[1]
                cb_1 = self.pdb_f.select('name CB and resindex ' + str(resind_1))
                ang1 = pr.calcAngle(ca_0, cb_0, cb_1)
                ang2 = pr.calcAngle(ca_1, cb_1, cb_0)
                if (ang1 < 80) or (ang2 < 80):
                    filter = False
            if filter:
                self.contacts_f.append((resnum_0, resnum_1))
        self.contacts_r = list()
        for nbr in self.nbrs_r:
            resnum_0 = nbr[0].getResnum()
            resnum_1 = nbr[1].getResnum()
            filter = True
            if np.abs(resnum_0 - resnum_1) > 6:
                resind_0 = nbr[0].getResindex()
                resind_1 = nbr[1].getResindex()
                ca_0 = nbr[0]
                cb_0 = self.pdb_r.select('name CB and resindex ' + str(resind_0))
                ca_1 = nbr[1]
                cb_1 = self.pdb_r.select('name CB and resindex ' + str(resind_1))
                ang1 = pr.calcAngle(ca_0, cb_0, cb_1)
                ang2 = pr.calcAngle(ca_1, cb_1, cb_0)
                if (ang1 < 80) or (ang2 < 80):
                    filter = False
            if filter:
                self.contacts_r.append((resnum_0, resnum_1))
            
    def initialize_sequence(self):
        """Sets the sequence to all neutral residues."""
        for i in self.surf_rns:
            if i in self.constrained_surf_rns_f:
                self.seq[i] = self.constrained_surf_rns_vals_f[i]
            else:
                self.seq[i] = 0

    def calc_gap(self):
        """Returns the energy gap of forward and reverse sequences."""
        e1 = np.sum(self.contact_En(con) for con in self.contacts_f)
        e2 = np.sum(self.contact_En(con) for con in self.contacts_r)
        return e1 - e2

    def contact_En(self, con):
        """Returns the contact energy between two residues in seq."""
        return En[self.seq[con[0]]][self.seq[con[1]]]

    def calc_En_f(self):
        """Returns the energy of the sequence mapped onto the pdb with forward topology."""
        return np.sum(self.contact_En(con) for con in self.contacts_f)

    def run_mc(self, num_iterations=1000000, kt_en_gap=1, kt_en_f=1, kt_seq_rep_len=0.4):
        """Runs Monte Carlo simulation that simultaneously optimizes:
        1. energy gap between topologies (maximizes),
        2. energy of forward topology sequence (minimizes),
        3. average length of repeat sequences (i.e. continuous strings
        such as 11111 in the sequence. minimizes)"""
        self.initialize_sequence()
        en_gap_old = 0
        en_f_old = 0
        seq_rep_len_old = len(self.surf_rns)
        for _ in range(num_iterations):
            key = random.choice(self.surf_rns)
            if key in self.constrained_surf_rns_f:  #accounts for constraints
                continue
            val = random.choice([0, 1, -1])
            old_val = self.seq[key]
            self.seq[key] = val
            en_gap_new = self.calc_gap()
            p_en_gap_new = np.exp((en_gap_old - en_gap_new) / kt_en_gap)
            en_f_new = self.calc_En_f()
            p_en_f_new = np.exp((en_f_old - en_f_new) / kt_en_f)
            seq_rep_len_new = np.mean([len(list(g)) for k, g in groupby(self.seq.values())])
            p_m_new = np.exp((seq_rep_len_old - seq_rep_len_new) / kt_seq_rep_len)
            if (p_en_f_new >= np.random.rand()) and (p_en_gap_new >= np.random.rand()) \
                    and (p_m_new >= np.random.rand()):
                en_gap_old = en_gap_new
                en_f_old = en_f_new
                seq_rep_len_old = seq_rep_len_new
                self.en_gaps.append(en_gap_new)
                self.en_fs.append(en_f_new)
                self.seq_rep_lens.append(seq_rep_len_new)
                self.seqs.append(self.seq.copy())
            else:
                self.seq[key] = old_val

    def find_pareto_front(self):
        """Returns the indices of the sequences that are in the pareto front,
        considering energy gap, energy of forward topology, and sequence repeat lengths."""
        costs = np.array(list(zip(self.en_gaps, self.en_fs, self.seq_rep_lens)))
        self.pareto_front = self.is_pareto_efficient_indexed(costs)

    @staticmethod
    def is_pareto_efficient_indexed(costs, return_mask=False):
        """
        :param costs: An (n_points, n_costs) array
        :param return_mask: True to return a mask, False to return integer indices of efficient points.
        :return: An array of indices of pareto-efficient points.
            If return_mask is True, this will be an (n_points, ) boolean array
            Otherwise it will be a (n_efficient_points, ) integer array of indices.

        This code is from username Peter at
        https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
        """
        is_efficient = np.arange(costs.shape[0])
        n_points = costs.shape[0]
        next_point_index = 0  # Next index in the is_efficient array to search for

        while next_point_index < len(costs):
            nondominated_point_mask = np.any(costs <= costs[next_point_index], axis=1)
            is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
            costs = costs[nondominated_point_mask]
            next_point_index = np.sum(nondominated_point_mask[:next_point_index]) + 1

        if return_mask:
            is_efficient_mask = np.zeros(n_points, dtype=bool)
            is_efficient_mask[is_efficient] = True
            return is_efficient_mask
        else:
            return is_efficient

    def find_nearest_utopian_pt(self, weight_en_gap=1, weight_en_f=1, weight_seq_rep_len=1):
        en_gaps_pareto = [self.en_gaps[i] for i in self.pareto_front]
        en_fs_pareto = [self.en_fs[i] for i in self.pareto_front]
        seq_rep_lens_pareto = [self.seq_rep_lens[i] for i in self.pareto_front]
        max_en_gap = max(en_gaps_pareto)
        max_en_f = max(en_fs_pareto)
        max_seq_rep_len = max(seq_rep_lens_pareto)
        min_en_gap = min(en_gaps_pareto)
        min_en_f = min(en_fs_pareto)
        min_seq_rep_len = min(seq_rep_lens_pareto)

        pareto_points = np.array([((self.en_gaps[i] - max_en_gap) / (min_en_gap - max_en_gap),
                                   (self.en_fs[i] - max_en_f) / (min_en_f - max_en_f),
                                   (self.seq_rep_lens[i] - max_seq_rep_len) / (min_seq_rep_len - max_seq_rep_len))
                                  for i in self.pareto_front])
        dists = cdist(pareto_points, np.array([[weight_en_gap, weight_en_f, weight_seq_rep_len]]))
        self.nearest_utopian_pt = self.pareto_front[(dists == min(dists)).flatten()]
        if len(self.nearest_utopian_pt) > 1:
            self.nearest_utopian_pt = self.nearest_utopian_pt[0]
        else:
            self.nearest_utopian_pt = int(self.nearest_utopian_pt)

    def map_seq_resnums_to_pdb(self, pdb):
        """Maps resnums of sequence to resnums of pdb."""
        self.resnum_conv = dict()
        self.resnum_conv_reverse = dict()
        for rn in self.seq.keys():
            sel = self.pdb_f.select('name CA and resnum ' + str(rn))
            pdb_rn = pdb.select('within 0.05 of sel', sel=sel).getResnums()[0]
            self.resnum_conv[rn] = pdb_rn
            self.resnum_conv_reverse[pdb_rn] = rn

    def set_charge_groups(self, seq):
        """Sets charged groups from sequence seq, which may be found via find_top(n)."""
        self.negs = list()
        self.neuts = list()
        self.poss = list()
        for rn, q in seq.items():
            if q == -1:
                self.negs.append(self.resnum_conv[rn])
            if q == 0:
                self.neuts.append(self.resnum_conv[rn])
            if q == 1:
                self.poss.append(self.resnum_conv[rn])

    def print_charge_groups(self):
        """Prints pymol resnum selection strings of the residues mapped to pdb
        via map_seq_resnums_to_pdb."""
        print('neutral residues= ' + '+'.join([str(i) for i in self.neuts]))
        print('negative residues= ' + '+'.join([str(i) for i in self.negs]))
        print('positive residues= ' + '+'.join([str(i) for i in self.poss]))

    def save_sequence(self, seq, outdir='./', filetag=''):
        if outdir[-1] != '/':
            outdir += '/'
        with open(outdir + 'surface_sequence' + filetag + '.txt', 'w') as outfile:
            outfile.write('resnum charge \n')
            for rn, q in seq.items():
                outfile.write(str(self.resnum_conv[rn]) + ' ' + str(q) + ' \n')


def set_bonds(prody_pdb):
    """Sets backbone bonds of chain based on proximity of atoms."""
    bb_sel = prody_pdb.select('protein and name N C CA')
    dm = pr.buildDistMatrix(bb_sel)
    ind = np.where((np.tril(dm) < 1.7) & (np.tril(dm) > 0))
    atom_ind = bb_sel.getIndices()
    prody_pdb.setBonds([(atom_ind[i], atom_ind[j]) for i, j in zip(ind[0], ind[1])])


def check_helix_lengths(hels):
    """Checks the helices that need to be equal length"""
    hels = list(hels)
    num_helices = len(hels)
    f = np.arange(1, num_helices)
    r = np.array(list(reversed(range(1, num_helices))))
    eq_len_helices = r != f
    print(eq_len_helices)
    lens = np.array([len(hel) for hel in hels[1:]])
    print(lens)
    if len(np.unique(lens[eq_len_helices])) == 1:
        return True
    else:
        return False


def order_helices(hels):
    """Cyclically permutes the helices until they achieve the right length"""
    hels = deque(hels)
    i = 0
    while not check_helix_lengths(hels):
        hels.rotate(1)
        i += 1
        if i > 10:
            break
            raise ValueError('Helices are not same length.')
    return list(hels)


def number_helices(pdb, reverse=True):
    """Renumbers a prody pdb/selection for forward or reverse topology.
    Returns a new prody object that is renumbered."""
    pdb = pdb.copy()
    set_bonds(pdb)
    i_end = 0
    hels = list()
    hel_inds = set()
    for i in pdb.iterBonds():
        if np.abs(i.getIndices()[0] - i_end) > 1:
            hels.append(hel_inds)
            hel_inds = set()
            i_end = i.getIndices()[-1]
            continue
        hel_inds |= set([b for b in i.getIndices()])
        i_end = i.getIndices()[-1]
    hels.append(hel_inds)
    hels = order_helices(hels)
    if reverse:
        order = [0]
        for n in list(reversed(range(1, len(hels)))):
            order.append(n)
    else:
        order = list(range(len(hels)))
    j = 1
    all_resnums = list()
    for o in order:
        resnums = list()
        hel = hels[o]
        ris = sorted(set(pdb.select('index ' + ' '.join([str(i) for i in hel])).getResindices()))
        for ri in ris:
            hel_sel = pdb.select('resindex ' + str(ri))
            resnums.extend(len(hel_sel) * [j])
            j += 1
        all_resnums.append(resnums)
    new_resnums = list()
    for o in order:
        new_resnums.extend(all_resnums[o])
    pdb.setResnums(new_resnums)
    return pdb

#to make dssp file
def parse_dssp(dssp_file):
    """parses a DSSP file and returns a dictionary of
    keys=(resnum, chain) and values=(dssp code)"""
    dssp_ss = dict()
    with open(dssp_file, 'r') as infile:
        for line in infile:
            if line.startswith('  #  RESIDUE'):
                break
        for line in infile:
            if line[13] == '!':
                continue
            resnum = int(line[5:10])
            chain = line[11]
            if line[16] != ' ':
                dssp_ss[(resnum, chain)] = line[16]
            else:
                dssp_ss[(resnum, chain)] = '-'
    return dssp_ss

from scipy.spatial import Delaunay
import prody as pr
import numpy as np
# import itertools
#from .pointTriangleDistance import pointTriangleDistance as distance
from numba import jit
import copy
# from scipy.optimize import linprog


# No numba version:
# def vol(a, b, c, d):
#     return np.abs(np.linalg.det(np.array([(a - d), (b - d), (c - d)]))) / 6

@jit("f8(f8[:],f8[:],f8[:],f8[:])", nopython=True, cache=True)
def vol(a, b, c, d):
    M = np.zeros((3, 3))
    M[0, :] = np.subtract(a, d)
    M[1, :] = np.subtract(b, d)
    M[2, :] = np.subtract(c, d)
    return np.abs(np.linalg.det(M)) / 6


@jit("f8(f8[:,:])", nopython=True, cache=True)
def get_radius(points):
    a = np.linalg.norm(points[0] - points[1])
    a1 = np.linalg.norm(points[2] - points[3])
    b = np.linalg.norm(points[0] - points[2])
    b1 = np.linalg.norm(points[1] - points[3])
    c = np.linalg.norm(points[0] - points[3])
    c1 = np.linalg.norm(points[1] - points[2])
    p = (a * a1 + b * b1 + c * c1) / 2
    V = vol(points[0], points[1], points[2], points[3])
    if V > 0:
        return 1 / (6 * V) * np.sqrt(p * (p - a * a1) * (p - b * b1) * (p - c * c1))
    else:
        return np.inf


@jit("i4[:,:](f8[:,:], i4[:,:], f8)", nopython=True, cache=True)
def _calc_alpha_simplex(C, S, a):
    M = S.shape[0]
    N = S.shape[1]
    Result = np.zeros((M, N))
    j = 0
    for i in range(M):
        s = S[i, :]
        ps = C[s]
        r = get_radius(ps)
        if r < a:
            Result[j, :] = s
            j += 1
    return Result[:j, :].astype(np.int32)


combos = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])

@jit("i4[:,:](i4[:,:], i8[:,:])", nopython=True, cache=True)
def make_simplex_set(S, combos):
    M = S.shape[0] * 4
    N = S.shape[1] - 1
    R = np.zeros((M, N), dtype=np.int32)
    for k, s in enumerate(range(0, M, 4)):
        for i in range(4):
            for j in range(3):
                R[s + i, j] = S[k, combos[i, j]]
    return R


def normal(ps):
    v1 = ps[1] - ps[0]
    v2 = ps[2] - ps[0]
    crossprod = np.cross(v1, v2)
    return crossprod / np.linalg.norm(crossprod)


class AlphaHull:

    def __init__(self, alpha=9):
        self.alpha = alpha
        self.hull = None
        self._hull = None
        self.tri = None
        self._tri = None
        self.coords = None
        self.simplices = None
        self.resindices = None
        # self.hull_points = None
#        self.hull_points_tr = None
#        self.len_hull_points = None

    def set_coords(self, pdb):
        type1 = isinstance(pdb, pr.atomic.selection.Selection)
        type2 = isinstance(pdb, pr.atomic.atomgroup.AtomGroup)
        if type1 or type2:
            self._set_coords(pdb)
        elif isinstance(pdb, np.ndarray):
            self.coords = pdb
        else:
            raise ValueError('*pdb* must be prody instance or numpy array')

    def _set_coords(self, pdb):
        """pdb is a prody object. pdb should have CB atoms where appropriate."""
        self.coords = pdb.select('name CB or (resname GLY and name CA)').getCoords()

    def set_tri(self):
        self.tri = Delaunay(self.coords)
        self._tri = copy.deepcopy(self.tri)

    def set_resindices(self, pdb):
        """pdb is a prody object. pdb should have CB atoms where appropriate."""
        self.resindices = pdb.select('name CB or (resname GLY and name CA)').getResindices()

    def calc_alpha_simplices(self):
        if self.tri is None:
            self.set_tri()

        self.tri.simplices.sort() #= np.sort(self.tri.simplices)
        self.tri.simplices = self.tri.simplices[self.tri.simplices[:, 0].argsort()]

        # numba compiled version is twice as fast:
        self.simplices = _calc_alpha_simplex(self.coords, self.tri.simplices, self.alpha)
        self._tri.simplices = self.simplices
        self._tri.neighbors = self.simplices
        #non compiled version:
        # self.simplices = np.array([simplex for simplex in self.tri.simplices
        #                            if get_radius(self.coords[simplex]) < self.alpha], dtype=np.int32)
        # self.simplices = np.array([simplex for simplex in sorted([sorted(tuple(s)) for s in self.tri.simplices])
        #                            if get_radius(self.coords[simplex]) < self.alpha], dtype=np.int32)

    def calc_hull(self):
        if self.simplices is None:
            self.calc_alpha_simplices()

        # simpl_set = [s[list(i)] for s in self.simplices
        #              for i in itertools.combinations(range(4), 3)]
        simpl_set = make_simplex_set(self.simplices, combos)

        un, ind, co = np.unique(simpl_set, axis=0,
                                return_counts=True, return_index=True)
        self.hull = np.array([simpl_set[i] for i in ind[co == 1]], dtype=np.int32)
        # self.hull_points = self.coords[list(set(self.hull.flatten()))]
        # self._hull = Delaunay(self.hull_points)
        # self.hull_points_tr = self.hull_points.T
        # self.len_hull_points = self.hull_points.shape[0]

    # def pnt_in_hull(self, pnt):
    #     '''
    #     Checks if `pnt` is inside the alpha hull.
    #     `pnt` -- point array of shape (3,)
    #     '''
    #
    #     if self.hull is None:
    #         self.calc_hull()
    #
    #     alph = AlphaHull(self.alpha)
    #     alph.set_coords(np.concatenate((self.coords, [pnt])))
    #     alph.calc_hull()
    #
    #     if np.array_equal(alph.hull, self.hull):
    #         return 1
    #     return -1

    # def pnt_in_hull(self, pnt):
    #     return self._tri.find_simplex(pnt) >= 0

       # if self.hull is None:
       #     self.calc_hull()
       #
       # c = np.zeros(self.len_hull_points)
       # A = np.r_[self.hull_points_tr, np.ones((1, self.len_hull_points))]
       # b = np.r_[pnt, np.ones(1)]
       # lp = linprog(c, A_eq=A, b_eq=b)
       # return lp.success

    def pnts_in_hull(self, pnts):
        return self._tri.find_simplex(pnts) >= 0

    # def pnts_in_hull_threshold(self, pnts, percent_buried):
    #     num_pts = len(pnts)
    #     in_out = list()
    #     num_out = 0
    #     for pnt in pnts:
    #         in_hull = self.pnt_in_hull(pnt)
    #         if not in_hull:
    #             num_out += 1
    #         if num_out/num_pts > (1 - percent_buried):
    #             return False, in_out
    #         in_out.append(in_hull)
    #     return True, in_out

    def get_pnt_distance(self, pnt):
        distances = []
        inout = self.pnts_in_hull(pnt)
        if inout:
            minmax = min
            inout = 1
        else:
            minmax = max
            inout = -1
        for i in range(len(self.hull)):
            distances.append(inout * distance(pnt, self.coords[self.hull[i]]))
        return minmax(distances)

    def get_pnts_distance(self, pnts):
        return [self.get_pnt_distance(pnt) for pnt in pnts]


def partition_res_by_burial(pdb_ala, alpha=9):
    """Returns residue indices of exposed, intermediate, and buried residues
    based on CA hull and CB hull."""
    ahull_ca = AlphaHull(alpha=alpha)
    ahull_ca.coords = pdb_ala.select('name CA').getCoords()
    ahull_ca.calc_hull()
    ahull_cb = AlphaHull(alpha=alpha)
    ahull_cb.set_coords(pdb_ala)
    ahull_cb.calc_hull()
    ahull_cb.set_resindices(pdb_ala)
    cb_in_ca_hull = ahull_ca.pnts_in_hull(ahull_cb.coords)
    resindices_cb_in_ca_hull = set(ahull_cb.resindices[cb_in_ca_hull])
    resindices_cb_hull = set(ahull_cb.resindices[np.unique(ahull_cb.hull)])
    resindices_not_cb_hull = set(ahull_cb.resindices) - resindices_cb_hull
    resindices_exposed = resindices_cb_hull - resindices_cb_in_ca_hull
    resindices_intermediate = resindices_cb_in_ca_hull - resindices_not_cb_hull
    resindices_buried = resindices_cb_in_ca_hull - resindices_intermediate
    res_ = resindices_not_cb_hull - resindices_buried
    resindices_intermediate |= res_
    return resindices_exposed, resindices_intermediate, resindices_buried

#
#
#
# def pnt_in_cvex_hull(hull, pnt):
#     '''
#     Checks if `pnt` is inside the convex hull.
#     `hull` -- a QHull ConvexHull object
#     `pnt` -- point array of shape (3,)
#     '''
#     new_hull = ConvexHull(np.concatenate((hull.points, [pnt])))
#     if np.array_equal(new_hull.vertices, hull.vertices):
#         return 1
#     return -1
#
# def pnt_in_alpha_hull(hull, pnt, alpha=7):
#     '''
#     Checks if `pnt` is inside the convex hull.
#     `hull` -- a QHull ConvexHull object
#     `pnt` -- point array of shape (3,)
#     '''
#     new_hull = get_alpha_hull(np.concatenate((hull.points, [pnt])), alpha=7)
#     if np.array_equal(new_hull.vertices, hull.vertices):
#         return 1
#     return -1
#
#
# def normal(ps):
#     v1 = ps[1] - ps[0]
#     v2 = ps[2] - ps[0]
#     crossprod = np.cross(v1, v2)
#     return crossprod / np.linalg.norm(crossprod)
#
#
# def distance(p, ps):
#     return abs(np.dot(normal(ps), p - ps[0]))
#
#
# def get_hull_coords(pdb):
#     """pdb is a prody object"""
#     coords = pdb.select('name CB or (resname GLY and name CA)').getCoords()
#     return coords
#
#
# def get_tri(coords):
#     return Delaunay(coords)
#
#
# def get_hull(simplices):
#     simpl_set = [s[list(i)] for s in simplices for i in itertools.combinations(range(4), 3)]
#     un, ind, co = np.unique(simpl_set, axis=0, return_counts=True, return_index=True)
#     hull = np.array([simpl_set[i] for i in ind[co == 1]], dtype=np.int32)
#     return hull
#
#
# def get_alpha_simplices(alpha, tri):
#     return np.array(list(filter(lambda simplex: get_radius(tri.points[simplex]) < alpha, tri.simplices)),
#                     dtype=np.int32)
#
#
#
#
#
# def get_alpha_hull(pdb, alpha=7):
#     type1 = isinstance(pdb, pr.atomic.selection.Selection)
#     type2 = isinstance(pdb, pr.atomic.atomgroup.AtomGroup)
#     if type1 or type2:
#         coords = get_hull_coords(pdb)
#     elif isinstance(pdb, np.ndarray):
#         coords = pdb
#
#     tri = get_tri(coords)
#     simplices = get_alpha_simplices(alpha, tri)
#     return get_hull(simplices)



#########################################################
# Example for getting  minimum distances of all CB atoms in a protein to the convex hull (protein surface)
# consisting of all CB coordinates of that protein.  This is the relevant function for the database distances.
#
# pdb_ala = pr.parsePDB('protein_all_ALA.pdb')
# cb_coords = pdb_ala.select('name CB and chain A').getCoords()    # Make sure these selections
# resnums = pdb_ala.select('name CB and chain A').getResnums()     # are only for the relevant chain.
# tri = Delaunay(cb_coords)
# min_distances = []
# for p in tri.points:
#     distances = [distance(p, tri.points[tri.convex_hull[i]])
#                  for i in range(len(tri.convex_hull))]
#     min_distances.append(min(distances))
# resnums_distances = list(zip(resnums, min_distances))


#######################################################
# Example for getting  minimum distance of a query point [x,y,z] to the convex hull (protein surface).
# The distance is negative if the point lies *outside* the hull.  (Positive if inside)
#
# hull = ConvexHull(cb_coords)
# pdb_w_ligand = pr.parsePDB('protein_with_ligand.pdb')
# min_distances_apx = []
# for p in pdb_w_ligand.select('resname APX').getCoords():
#     distances = []
#     inout = pnt_in_cvex_hull(hull, p)
#     if inout == 1:
#         minmax = min
#     else:
#         minmax = max
#     for i in range(len(tri.convex_hull)):
#         distances.append(inout * distance(p, tri.points[tri.convex_hull[i]]))
#     min_distances_apx.append(minmax(distances))
# names_distances = list(zip(pdb_w_ligand.select('resname APX').getNames(), min_distances_apx))


# This was taken from joshuashaffer @ https://gist.github.com/joshuashaffer/99d58e4ccbd37ca5d96e
# with minimal changes

#!/usr/bin/env python
#
# Tests distance between point and triangle in 3D. Aligns and uses 2D technique.
#
# Was originally some code on mathworks

# import numpy
from numpy import dot
from math import sqrt



def pointTriangleDistance(P, TRI):
    # function [dist,PP0] = pointTriangleDistance(TRI,P)
    # calculate distance between a point and a triangle in 3D
    # SYNTAX
    #   dist = pointTriangleDistance(TRI,P)
    #   [dist,PP0] = pointTriangleDistance(TRI,P)
    #
    # DESCRIPTION
    #   Calculate the distance of a given point P from a triangle TRI.
    #   Point P is a row vector of the form 1x3. The triangle is a matrix
    #   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
    #   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
    #   to the triangle TRI.
    #   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
    #   closest point PP0 to P on the triangle TRI.
    #
    # Author: Gwolyn Fischer
    # Release: 1.0
    # Release date: 09/02/02
    # Release: 1.1 Fixed Bug because of normalization
    # Release: 1.2 Fixed Bug because of typo in region 5 20101013
    # Release: 1.3 Fixed Bug because of typo in region 2 20101014

    # Possible extention could be a version tailored not to return the distance
    # and additionally the closest point, but instead return only the closest
    # point. Could lead to a small speed gain.

    # Example:
    # %% The Problem
    # P0 = [0.5 -0.3 0.5]
    #
    # P1 = [0 -1 0]
    # P2 = [1  0 0]
    # P3 = [0  0 0]
    #
    # vertices = [P1; P2; P3]
    # faces = [1 2 3]
    #
    # %% The Engine
    # [dist,PP0] = pointTriangleDistance([P1;P2;P3],P0)
    #
    # %% Visualization
    # [x,y,z] = sphere(20)
    # x = dist*x+P0(1)
    # y = dist*y+P0(2)
    # z = dist*z+P0(3)
    #
    # figure
    # hold all
    # patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8)
    # plot3(P0(1),P0(2),P0(3),'b*')
    # plot3(PP0(1),PP0(2),PP0(3),'*g')
    # surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
    # view(3)

    # The algorithm is based on
    # "David Eberly, 'Distance Between Point and Triangle in 3D',
    # Geometric Tools, LLC, (1999)"
    # http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    #
    #        ^t
    #  \     |
    #   \reg2|
    #    \   |
    #     \  |
    #      \ |
    #       \|
    #        *P2
    #        |\
    #        | \
    #  reg3  |  \ reg1
    #        |   \
    #        |reg0\
    #        |     \
    #        |      \ P1
    # -------*-------*------->s
    #        |P0      \
    #  reg4  | reg5    \ reg6
    # rewrite triangle in normal form
    B = TRI[0, :]
    E0 = TRI[1, :] - B
    # E0 = E0/sqrt(sum(E0.^2)); %normalize vector
    E1 = TRI[2, :] - B
    # E1 = E1/sqrt(sum(E1.^2)); %normalize vector
    D = B - P
    a = dot(E0, E0)
    b = dot(E0, E1)
    c = dot(E1, E1)
    d = dot(E0, D)
    e = dot(E1, D)
    f = dot(D, D)

    #print "{0} {1} {2} ".format(B,E1,E0)
    det = a * c - b * b
    s = b * e - c * d
    t = b * d - a * e

    # Terible tree of conditionals to determine in which region of the diagram
    # shown above the projection of the point into the triangle-plane lies.
    if (s + t) <= det:
        if s < 0.0:
            if t < 0.0:
                # region4
                if d < 0:
                    t = 0.0
                    if -d >= a:
                        s = 1.0
                        sqrdistance = a + 2.0 * d + f
                    else:
                        s = -d / a
                        sqrdistance = d * s + f
                else:
                    s = 0.0
                    if e >= 0.0:
                        t = 0.0
                        sqrdistance = f
                    else:
                        if -e >= c:
                            t = 1.0
                            sqrdistance = c + 2.0 * e + f
                        else:
                            t = -e / c
                            sqrdistance = e * t + f

                            # of region 4
            else:
                # region 3
                s = 0
                if e >= 0:
                    t = 0
                    sqrdistance = f
                else:
                    if -e >= c:
                        t = 1
                        sqrdistance = c + 2.0 * e + f
                    else:
                        t = -e / c
                        sqrdistance = e * t + f
                        # of region 3
        else:
            if t < 0:
                # region 5
                t = 0
                if d >= 0:
                    s = 0
                    sqrdistance = f
                else:
                    if -d >= a:
                        s = 1
                        sqrdistance = a + 2.0 * d + f;  # GF 20101013 fixed typo d*s ->2*d
                    else:
                        s = -d / a
                        sqrdistance = d * s + f
            else:
                # region 0
                invDet = 1.0 / det
                s = s * invDet
                t = t * invDet
                sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f
    else:
        if s < 0.0:
            # region 2
            tmp0 = b + d
            tmp1 = c + e
            if tmp1 > tmp0:  # minimum on edge s+t=1
                numer = tmp1 - tmp0
                denom = a - 2.0 * b + c
                if numer >= denom:
                    s = 1.0
                    t = 0.0
                    sqrdistance = a + 2.0 * d + f;  # GF 20101014 fixed typo 2*b -> 2*d
                else:
                    s = numer / denom
                    t = 1 - s
                    sqrdistance = s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f

            else:  # minimum on edge s=0
                s = 0.0
                if tmp1 <= 0.0:
                    t = 1
                    sqrdistance = c + 2.0 * e + f
                else:
                    if e >= 0.0:
                        t = 0.0
                        sqrdistance = f
                    else:
                        t = -e / c
                        sqrdistance = e * t + f
                        # of region 2
        else:
            if t < 0.0:
                # region6
                tmp0 = b + e
                tmp1 = a + d
                if tmp1 > tmp0:
                    numer = tmp1 - tmp0
                    denom = a - 2.0 * b + c
                    if numer >= denom:
                        t = 1.0
                        s = 0
                        sqrdistance = c + 2.0 * e + f
                    else:
                        t = numer / denom
                        s = 1 - t
                        sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f

                else:
                    t = 0.0
                    if tmp1 <= 0.0:
                        s = 1
                        sqrdistance = a + 2.0 * d + f
                    else:
                        if d >= 0.0:
                            s = 0.0
                            sqrdistance = f
                        else:
                            s = -d / a
                            sqrdistance = d * s + f
            else:
                # region 1
                numer = c + e - b - d
                if numer <= 0:
                    s = 0.0
                    t = 1.0
                    sqrdistance = c + 2.0 * e + f
                else:
                    denom = a - 2.0 * b + c
                    if numer >= denom:
                        s = 1.0
                        t = 0.0
                        sqrdistance = a + 2.0 * d + f
                    else:
                        s = numer / denom
                        t = 1 - s
                        sqrdistance = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f

    # account for numerical round-off error
    if sqrdistance < 0:
        sqrdistance = 0

    dist = sqrt(sqrdistance)

    #PP0 = B + s * E0 + t * E1
    return dist  #, PP0

# if __name__ == '__main__':
#     TRI = numpy.array([[0.0,-1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,0.0]])
#     P = numpy.array([0.5,-0.3,0.5])
#     dist, pp0 = pointTriangleDistance(TRI,P)
#     print dist
#     print pp0