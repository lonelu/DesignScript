import prody as pr
import os

file_path = '/mnt/e/GitHub_Design/DesignScript/data/sasa/'

file1 = 'full_7a_10a.pdb'
pdb_prody1 = pr.parsePDB(file_path + file1)
ca1 = pdb_prody1.select('name CA')

#ca_coords = pdb_prody.select('name CA').getCoords()

file2 = 'full_7a_14a.pdb'
pdb_prody2 = pr.parsePDB(file_path + file2)
ca2 = pdb_prody2.select('name CA')


import matplotlib.pyplot as plt


pr.showContactMap(ca1, cmap = 'Reds')
plt.savefig(file_path + os.path.basename(file1) + '.png')
plt.close()

pr.showContactMap(ca2, cmap = 'Reds')
plt.savefig(file_path + os.path.basename(file2) + '.png')
plt.close()