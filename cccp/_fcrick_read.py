import os
import sys


def get_all_par_file(workdir):
    all_par_files = []

    for par_path in os.listdir(workdir):
        if not par_path.endswith(".par"):
            continue
        all_par_files.append(workdir + par_path)
    return all_par_files

def read_par_info(par_path):
    thedict = {'name':par_path.split('/')[-1]}
    with open(par_path, 'r') as f:
        for r in f.readlines():
            ss = r.split(' = ')
            if ':' in ss[0]: 
                key = ss[0].split(': ')[1]
            else:
                key = ss[0]
            value = ss[1].split('\n')[0]
            if ' ' in value:
                value = value.split(' ')[0]
            thedict[key] = value
    return thedict

def write_par_summary(workdir, outdir, out_name = '_summary_par_info.txt'):
    '''
    read all par files, write to a tab deliminated file.
    '''
    par_dicts = []
    files = get_all_par_file(workdir)
    for f in files:
        _dict = read_par_info(f)
        par_dicts.append(_dict)

    with open(outdir + out_name, 'w') as f:
        f.write('\t'.join([k for k in list(par_dicts[0].keys())]) + '\n')
        for d in par_dicts:
            f.write('\t'.join([k for k in list(d.values())]) + '\n')
    


### Get single par info.
#par_path = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2_coil/_z_fix_c2/'
#thedict = read_par_info(par_path)

### Get all par info
def main():
    workdir = '/mnt/e/DesignData/ligands/LigandBB/_reverse_design/c2_coil/'
    indir = workdir + '_z_fix_c2/'
    outdir = workdir + '_z_fix_c2_enriched/'
    os.makedirs(outdir, exist_ok=True)

    write_par_summary(indir, outdir)
    return 

if __name__ == "__main__":
    main()

