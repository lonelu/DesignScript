pdbs = get_all_pbd_prody(workdir + '2_seq_cores_reps/')

pdb_prody = [p for p in pdbs if p.getTitle()=='3rqt_1'][0] 

nis = pdb_prody.select(metal_sel)

ni = nis[0]

all_aa = pdb_prody.select(aa + ' and within 2.6 of index ' + str(ni.getIndex()))

inds = np.unique(all_aa.getResindices())

ind = inds[1]

ext_inds = extend_res_indices([ind], pdb_prody, extend =1)

pdb_prody.getResindices()