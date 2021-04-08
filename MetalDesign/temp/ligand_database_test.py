
# Superimpose metal (for test purpose)
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

# get_metal_core_seq (for test purpose)
'''
import os

import prody as pr

metal_sel = 'name NI'

workdir = "/mnt/e/DesignData/ligands/NI/pfam2pdbs_sel/"

file_path = workdir + '1a5n4.pfam.pdb'



try:
    pdb_prody = pr.parsePDB(file_path)
except:
    print('not sure')

nis = pdb_prody.select(metal_sel)

metal_cores = []
count = 0

#for ni in nis:
ni = nis[0]

ni_index = ni.getIndex()
all_near = pdb_prody.select('nitrogen or oxygen').select('not water and within 2.6 of index ' + str(ni_index))
#if len(all_near) < 3:
#    continue    
inds = all_near.getResindices()
all_near_res = pdb_prody.select('name CA and resindex ' + ' '.join([str(ind) for ind in inds]))
#if len(all_near_res) < 1:
#    continue     
ext_inds = extend_res_indices(all_near_res.getResindices(), pdb_prody)
count += 1
sel_pdb_prody = pdb_prody.select('resindex ' + ' '.join([str(ind) for ind in ext_inds]) + 'and index '+ str(ni_index))
metal_cores.append((os.path.basename(file_path) + '_' + str(count), sel_pdb_prody))
'''

# test extend_res_indices
"""
inds_near_res = all_near_res.getResindices()

extend_inds = []
inds = set()
for ind in inds_near_res:
    for i in range(-2, 3):
        the_ind = ind + i
        if the_ind not in inds and the_ind>= 0 and connectivity_filter(pdb_prody, ind, the_ind):
            extend_inds.append(the_ind)
            inds.add(the_ind)
"""


#Test extract_all_core_seq
'''
import os
import prody as pr
import ligand_database as ld

metal_sel = 'name NI'

workdir = "/mnt/e/DesignData/ligands/NI/pfam2pdbs_sel/"

file_path = workdir + '1a5n4.pfam.pdb'

metal_cores = ld.get_metal_core_seq(file_path, metal_sel)

sub_workdir = 'seq_cores/'

cores = []
first_get = False
first = None

core = metal_cores
    
if not first_get:
    first = core[0]
    first_get = True
first[1].getIndices() 

if not os.path.exists(workdir + sub_workdir):
    os.mkdir(workdir + sub_workdir)

for c in core:
    a_coords = first[1].select(metal_sel)[0].getCoords().reshape(1,3)
    b_coords = c[1].select(metal_sel)[0].getCoords().reshape(1,3)

    pr.calcTransformation(b_coords, a_coords).apply(c[1])

    outfile = c[0]
    pr.writePDB(workdir + sub_workdir + outfile + '_.pdb', c[1])

'''

# test reduce_dup
'''
workdir = "/mnt/e/DesignData/ligands/NI_rcsb/seq_cores/"

cores = get_all_pbd_prody(workdir)

pdbs = [cores[6], cores[7]]

i = 0
j = 1

dic = dict()
clusters = []
clustered = set()

dic[i] = pdbs[i].select('name CA').getCoords()
cluster = [i]
clustered.add(i)

len(pdbs[i]) != len(pdbs[j])

pr.calcRMSD(pdbs[i], pdbs[j])
pr.calcTransformation(pdbs[j].select('name CA'), pdbs[i].select('name CA')).apply(pdbs[j])
pr.calcRMSD(pdbs[i], pdbs[j])

pdbs = cores
metal_sel = 'name NI'
clusters = reduce_dup(pdbs, metal_sel)
'''