#To find the ligand, search it in the RCSB PDB database. 
#Then in the Search Results, select Tabular Report -> Summary Report -> Structure
#Custorm report 'Return to Create Custom Report.'

import os
import csv

#"Entry ID	Chain ID	Ligand ID	Ligand Formula	Ligand MW	Ligand Name	InChI	InChI Key	Ligand SMILES	Ligand of Interest"

all_lines = []
path = "/mnt/e/DesignData/ligands/NI_rcsb/"
for file in os.listdir(path):
    if file.endswith(".csv"):
        with open(path + file, 'r') as f:
            reader = csv.reader(f)
            for r in reader:
                all_lines.append(r)
with open(path + 'all.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(all_lines[0])
    for r in all_lines:
        if not r[0] == 'Entry ID':
            writer.writerow(r)


import os
import csv
from prody import *
path = "/mnt/e/DesignData/ligands/NI_rcsb/"
all_pdbs = []
with open(path + 'all_NI.csv', 'r') as f:
    reader = csv.reader(f) 
    title = next(reader)
    for r in reader:
        if (r[0] !='') and (',' not in r[4]) and (r[4]!= '') and (float(r[4]) <= 2.5):
            all_pdbs.append(r[0])

fetchPDBviaFTP(all_pdbs[0], compressed = False)
pathPDBFolder(path)
for p in all_pdbs:
    fetchPDBviaFTP(p, compressed = False) 
#Then unzip them in linux with:  
cd "/mnt/e/DesignData/ligands/NI_rcsb/"
gunzip *.gz

    
#The Bio.PDB format is weird.
import sys
from Bio.PDB import PDBList
download_pdb = PDBList()
download_pdb.retrieve_pdb_file(pdb_code = '4UWX', obsolete = False, pdir = path,  file_format = 'pdb') 
download_pdb.retrieve_pdb_file(pdb_codes = all_pdbs, obsolete = False, pdir = path,  file_format = 'pdb') 
    