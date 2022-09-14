'''
OmegaFold is very fast, which could be used to substitute Ab initial fold (robetta).
The script here is used to superimpose the OmegaFold output pdb to the rosetta pdb and obtain the rmsd.
'''

import os
import prody as pr

'''

'''

workdir = '/mnt/e/DesignData/Chemodrugs/HB_RUC_2nd_599/'

bb_prep = pr.parsePDB(workdir + 'bb_prep_599.pdb')

rosetta_indir = workdir + 'output_sel/'
rosetta_outdir = workdir + 'output_sel/rosetta_output_ana/'
os.makedirs(rosetta_outdir, exist_ok=True)

omega_indir = workdir + 'output_sel/omegafold_output/'
omega_outdir = workdir + 'output_sel/omegafold_output_ana/'
os.makedirs(omega_outdir, exist_ok=True)

rmsds = []

for _pdb_file in os.listdir(rosetta_indir):
    if not '.pdb' in _pdb_file:
        continue
    _pdb = pr.parsePDB(rosetta_indir + _pdb_file)
    pr.calcTransformation(_pdb.select('bb'), bb_prep.select('bb')).apply(_pdb)
    _rmsd = pr.calcRMSD(_pdb.select('bb'), bb_prep.select('bb'))

    _omgpdb = pr.parsePDB(omega_indir + _pdb_file.split('.')[0] + '_A.pdb')
    pr.calcTransformation(_omgpdb.select('bb'), _pdb.select('bb')).apply(_omgpdb)
    _rmsd_dl = pr.calcRMSD(_omgpdb.select('bb'), _pdb.select('bb'))

    pr.writePDB(rosetta_outdir + _pdb.getTitle() + '.pdb', _pdb)
    pr.writePDB(omega_outdir + _omgpdb.getTitle() + '.pdb', _omgpdb)
    rmsds.append((_pdb.getTitle(), _rmsd, _rmsd_dl))

with open(rosetta_indir + '_rmsd_summary.tsv', 'w') as f:
    f.write('Title\trmsd_Design_Rosetta\trmsd_Rosetta_DL\n')
    for r in rmsds:
        f.write('\t'.join([str(_r) for _r in r]) + '\n')


