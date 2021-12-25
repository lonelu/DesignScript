import numpy as np

filenames = ['helix_query.pdb', 'TYR_ASP_483_C=2.957_103139_1.pdb']

with open('helix_103139_1.pdb', 'w') as outfile:
    for names in filenames:
        with open(names) as infile:
            outfile.write(infile.read())
        outfile.write("\n")
