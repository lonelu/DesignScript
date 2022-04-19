'''
The script tried to align two pdbs.
'''
import os
import prody as pr

workdir = '/mnt/e/DesignData/_temp/yang/'
outdir = workdir + 'prody_out/'
os.makedirs(outdir)

pdbs = []
pdb_titles = []
for p in os.listdir(workdir):
    if '.pdb' not in p:
        continue
    pdb_titles.append(p)
    pdbs.append(pr.parsePDB(workdir + p))


rmsds = []
for i in range(0, len(pdbs)-1, 2):
    target = pdbs[i]
    mobile = pdbs[i+ 1]
    target_sel = target.select('serial 2 4 21 43 45 62 84 86 103')
    mobile_sel = mobile.select('serial 2 4 21 43 45 62 84 86 103')

    rmsd_bf = pr.calcRMSD(mobile, target)
    pr.calcTransformation(mobile_sel, target_sel).apply(mobile)
    pr.writePDB(outdir + pdbs[i+1].getTitle() + '_tr.pdb', pdbs[i+ 1])

    rmsd = pr.calcRMSD(mobile, target)
    rmsds.append((rmsd_bf, rmsd))