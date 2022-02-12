# annotate PDBs in final database from 

import sys, os
from collections import Counter

path2_nonRedundant_list 	= sys.argv[1]
path2headers 				= sys.argv[2]
path2expmethod_list 		= sys.argv[3]



## build up list of ~300 good pdbs
headerPathes	  = {}
nonRedundant_list = []
for i in open(path2_nonRedundant_list):
	if len(i) < 3: continue
	nonRedundant_list.append(i.rstrip())


expmethod 		 = []
expmethod_list	 = {}
for i in open(path2expmethod_list):
	if len(i) < 3: continue
	pdb, data = i[:4], ' '.join(i.split(' ')[1:-1])
	if pdb in nonRedundant_list:
		expmethod_list[pdb] = data
		data2 = ' '.join( data.rstrip().rsplit(' ')[1:-1] ) 
		expmethod.append(data2)
		print pdb, data
		headerPathes[pdb] = os.path.join( path2headers ,  '%s.txt' % pdb.upper() )

cnts = Counter(expmethod)
for k, v in cnts.items():
	print '#\t',k, v 


