from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.rdmolfiles import MolToPDBFile
from rdkit import rdBase
print(rdBase.rdkitVersion)

import os
import prody as pr

import py3Dmol
def drawit(m,p=None,confId=-1):
        mb = Chem.MolToMolBlock(m,confId=confId)
        if p is None:
            p = py3Dmol.view(width=400,height=400)
        p.removeAllModels()
        p.addModel(mb,'sdf')
        p.setStyle({'stick':{}})
        p.setBackgroundColor('0xeeeeee')
        p.zoomTo()
        return p.show()

workdir = '/mnt/e/DesignData/ligands_metal_prot/_lig_fe/'

#Load ligand method1.
m_eno = Chem.MolFromPDBFile(workdir + 'eno_fe2.pdb')
#Add Hs
m2=Chem.AddHs(m_eno)
AllChem.EmbedMolecule(m2)

#Load ligand method2.
m = Chem.MolFromSmiles('c1cc(ccc1CC(=O)C(=O)O)O')
m2=Chem.AddHs(m)
AllChem.EmbedMolecule(m2)

# generate confomers.
cids = AllChem.EmbedMultipleConfs(m2, numConfs=64)

rmslist = []
AllChem.AlignMolConformers(m2, RMSlist=rmslist)

from rdkit.Chem.rdmolfiles import MolToPDBFile
for i in range(64):
    MolToPDBFile(m2, workdir + 'all_ligs_rdkit_i3pa_test/i3pa_' + str(i) + '.pdb', confId = i)
