'''
After rosetta design, thousands of proteins scores are in seperate files. 
This script tried to combine them into one file for easier comparison. 
'''


'''
conda activate env_conda
python /mnt/e/GitHub_Design/DesignScript/Rosetta/combine_score_sc.py
'''

import os
import pandas
import shutil

_dir = '/wynton/home/degradolab/lonelu/DesignData/fluo/Rosetta/'
_dir = '/wynton/home/degradolab/lonelu/DesignData/Chemodrugs/NIR_2nd/'
workdir = _dir + 'output/'
outdir = _dir + 'output_sel/'


os.makedirs(outdir, exist_ok=True)

title = ''

all_lines = []
for file in os.listdir(workdir):
    if '.sc' not in file:
        continue
    with open(workdir + file, 'r') as f:
        for line in f.readlines():
            if 'SEQUENCE' in line:
                continue
            if 'total_score' in line and len(title) == 0:
                title = 'title'+ '\t' + '\t'.join([_l for _l in line.split(' ') if _l != ''])
            if 'total_score' in line:
                continue
            all_lines.append(file + '\t' + '\t'.join([_l for _l in line.split(' ') if _l != '']))

with open(outdir + '_summary.tsv', 'w') as f:
    f.write(title)
    for line in all_lines:
        f.write(line)


#>>> Copy the top 50 and align

pdb_file_dict = {}
for file in os.listdir(workdir):
    if not '.pdb' in file:
        continue
    id = file.split('.')[0].split('_')[-1]
    pdb_file_dict[id] = file

pd = pandas.read_csv(outdir + '_summary.tsv', sep='\t')

for title in pd.sort_values(by=['total_score'])['title'].iloc[0:30]:
    id = title.split('.')[0].split('_')[1]
    file = pdb_file_dict[id]
    shutil.copy(workdir + file, outdir + file)