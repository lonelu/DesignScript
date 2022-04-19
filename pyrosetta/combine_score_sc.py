import os

'''
python /mnt/e/GitHub_Design/DesignScript/pyrosetta/combine_score_sc.py
'''

#workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta_86-88-101/_rosetta_r1/output_1/'
workdir = '/mnt/e/DesignData/ligands/LigandBB/_lig_fe/_ntf2_rosetta_86-88-101/_rosetta_r3/output_59Y/'

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
                title = line
            if 'total_score' in line:
                continue
            all_lines.append(file + ' ' + line)

with open(workdir + '_summary.txt', 'w') as f:
    f.write('title ' + title)
    for line in all_lines:
        f.write(line)

