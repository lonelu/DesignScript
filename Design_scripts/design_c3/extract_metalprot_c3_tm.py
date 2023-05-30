'''
For Kehan's 6 helix bundles A-B-C-D-E-F, we want to find C3 binding from A-C-E or B-D-F. 
The script here is to extract the  A-C-E or B-D-F from all combinations in wynton.
'''


import os
import shutil

indir = '/wynton/home/degradolab/lonelu/DesignData/tm/queries/'
outdir = '/wynton/home/degradolab/lonelu/DesignData/tm/queries_extract/'

for dir_0 in os.listdir(indir):
    #print(indir + dir_0)
    if not os.path.isdir(indir + dir_0):
        continue
    for dir_1 in os.listdir(indir + dir_0):
        #print(indir + dir_0 + '/' + dir_1)
        if not os.path.isdir(indir + dir_0 + '/' + dir_1):
            #print('--------------------------')
            continue
        if dir_1 == 'represents':
            for _file_rep in os.listdir(indir + dir_0 + '/' + dir_1):
                xss = _file_rep.split('_')
                xs = xss[1].split('-')
                #print(xs)
                if (xs[0] == 'A' and xs[2]=='C' and xs[4]=='E') or (xs[0] == 'B' and xs[2]=='D' and xs[4]=='F'):
                    os.makedirs(outdir + dir_0 + '/' + dir_1 + '/', exist_ok=True)
                    shutil.copy(indir + dir_0 + '/' + dir_1 + '/' + _file_rep, outdir + dir_0 + '/' + dir_1 + '/' + _file_rep)
        else:
            xss = dir_1.split('_')
            xs = xss[1].split('-')
            #print(xs)
            if (xs[0] == 'A' and xs[2]=='C' and xs[4]=='E') or (xs[0] == 'B' and xs[2]=='D' and xs[4]=='F'):
                os.makedirs(outdir + dir_0 + '/' + dir_1 + '/', exist_ok=True)
                shutil.copytree(indir + dir_0 + '/' + dir_1, outdir + dir_0 + '/' + dir_1, dirs_exist_ok=True)