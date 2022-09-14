'''
After rosetta design, thousands of proteins scores are in seperate files. 
This script tried to combine them into one file for easier comparison. 
'''


'''
python /mnt/e/GitHub_Design/DesignScript/Rosetta/combine_score_sc.py
'''

#workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta_86-88-101/_rosetta_r1/output_1/'
#workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta_86-88-101/_rosetta_r3/output_59Y/'
#workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/design_0_rosetta821/output_1/'
#workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/design_1_rosetta/output_0/'
#workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig_318/loop_1_rosetta/output_0/'
#workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_10-14-139_HDH/Lig-0/design_0_rosetta/output_0/'
#workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_14-135-139_HHE/helix6_loop1_rosetta33/output_1/'
#workdir = '/mnt/e/DesignData/Metalloenzyme/HelixFe_TTS/helix6_14-135-139_HHE/helix6_loop2_rosetta/output_2/'
#workdir = '/mnt/e/DesignData/Metalloenzyme/1ukr/68-77-79_DHH/lig20/rosetta_design_0/output_0/'

import os
import pandas
import shutil

_dir = '/mnt/e/DesignData/Metalloenzyme/SAHA_Vorinostat/run_design_cgs3/SAHA_Rosetta20220729/rosetta/'
_dir = '/wynton/home/degradolab/lonelu/DesignData/Chemodrugs/HB_RUC_2nd_599/'
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


#>>> Copy the top 10 and align



pdb_file_dict = {}
for file in os.listdir(workdir):
    if not '.pdb' in file:
        continue
    id = file.split('.')[0].split('_')[-1]
    pdb_file_dict[id] = file

pd = pandas.read_csv(outdir + '_summary.tsv', sep='\t')

for title in pd.sort_values(by=['total_score'])['title'].iloc[0:10]:
    id = title.split('.')[0].split('_')[1]
    file = pdb_file_dict[id]
    shutil.copy(workdir + file, outdir + file)