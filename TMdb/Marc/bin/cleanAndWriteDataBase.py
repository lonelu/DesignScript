import os, sys, numpy as np, cPickle as pkl
from collections import defaultdict, Counter
from prody import * 

from tm_pdb_utilities import *

cPkl_file_path 	 	 = sys.argv[1]

master_db 			 = '40%%master-nr-TMdatabase_bychain/'
wholePDB_dir		 = '40%%IDfull_nr-TMpdbs/'
fasta_dir 			 = sys.argv[4]
OPM_dirpath			 = sys.argv[3]
all_Sequences 		 = []

if not os.path.exists( master_db ):
	os.mkdir( master_db )

if not os.path.exists( wholePDB_dir ):
	os.mkdir( wholePDB_dir )	

# load data, a directory of Pdb objects containing chain objects
dict_of_pdb_objs = pkl.load( open( cPkl_file_path , 'rb') )

# parse list for non-redundant pdb files
# during, print out 1) whole TM PDB 
# 2) each chain group representative & its interchain interactions
#
for i in open( sys.argv[2] ):
	pdb = dict_of_pdb_objs[ i.rstrip() ]

	print pdb, pdb.chainGroups.keys(), pdb.chainGroups.values()

	wholePDB_Path 	= os.path.join( wholePDB_dir, 'tm_%s.pdb' % pdb.pdb )
	all_chain_paths = [ os.path.join( OPM_dirpath, 'tm_%s.pdb' % x ) for x in sorted( pdb.chain_List.values() ) ]

	interChain_interactionsFound = []

	wholePDB_obj = 0
	for ch_path in all_chain_paths:

		chain = parsePDB( ch_path )
		# parse and write full pdb objects
		if not wholePDB_obj:
			wholePDB_obj = chain
		else:
			wholePDB_obj += chain

	writePDB( wholePDB_Path,  wholePDB_obj )

	# Now the whole PDB is loaded into memory, 
	#  Cycle through chain group representatives
	#  and lcoatte the novel inter-chain interactions
	for chain in pdb.chainGroups.keys():
		
		chain_rep = wholePDB_obj.select( 'chain %s ' % chain )

		redundant_chains = [ wholePDB_obj.select( 'chain %s ' % chX ) for chX in pdb.chainGroups[chain] if chain !=chX ]

		# locate interchain interactions by nearby
		nearby_residues = wholePDB_obj.select( 'same residue as (bb and exwithin 12 of chain %s)' % (chain ) )


		if nearby_residues:
			chain_rep += nearby_residues

		chain_rep_path = os.path.join( master_db, 'tm_%s_%s.pdb' % (pdb, chain ) )
	
		writePDB( chain_rep_path, chain_rep )


		### Keep track of the TM sequence counts for all amino acids
		# this info is in the fasta file for each 

		chain_fasta_path = os.path.join( fasta_dir, 'tm_%s_%s.fasta' % (pdb, chain )  )

		for seq in open( chain_fasta_path ):
			if seq[0] != '>':
				all_Sequences.extend( seq.rstrip () )


		

amino_acid_counts = Counter( all_Sequences )
total = np.sum( [amino_acid_counts.values()] )

for k, v in amino_acid_counts.items():
	print k, v, round( float( v) / total , 3 )



