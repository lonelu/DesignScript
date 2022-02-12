import os, sys, numpy as np
from collections import defaultdict
from itertools import combinations, product 
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import cPickle as pkl

### MAKE SURE tm_pdb_utilities.py is in your $PYTHON_PATH
from tm_pdb_utilities import *

subsMat = matlist.blosum62

# Look at each helical TM pdb, compare what chains are unique 
# (>90% of any TM segment, and >70% overall average pairwise idenitity)
# Choose the largest chain in the chain group 

# input (1) list good redundant set of PDBs, sorted by chain with % helical data
# input (2) path to directory for each TM chain's fasta file
# input (3) path to directory for each TM chain's pdb file

# output(1) list of all chain groups in each pdb, with noted representitive chains
# output(2 & 3) cPickle and JSON outputs holding the PDB objects

# example command line
# python bin/findChainGroups.py suitable_redundant_pdbByChain.txt selectedOPM/ selectedfasta/
# WRITTEN FOR PYTHON 2.7

### ...seriously ... MAKE SURE tm_pdb_utilities.py is in your $PYTHON_PATH





inputList_byChain 	= sys.argv[1]
pdb_dir				= sys.argv[2]
fasta_dir			= sys.argv[3]


dict_of_pdb_objs 		= {}
helical_fraction_data 	= []
average_pairwise_seqID  = []
chain_group_log_txt		= ''

beta_or_small_BADpdbs 	= []

# parse input file, 
for fin in open( inputList_byChain ):
	error_message ='\ninconsistant formatting of input, is pdb mentioned before pdb_chain?... %s quitting\n\n' % fin


	line = fin.split() 

	if len( line[0] ) ==  4:   		# handle pdb

		exp_type, resolution = ' '.join( line[1:-1] ), line[-1]
		dict_of_pdb_objs[ line[0] ] = Pdb( line[0] , exp_type, resolution , pdb_dir)

	elif len( line[0] ) == 6:		# 	handle chain

		# Skip non-helical segments
		helical_fraction = float( line[-1] )
		helical_fraction_data.append( helical_fraction )

		tm_length = float( line[1].split('=')[-1] )

		if helical_fraction < 60:
			continue

		try:	
			pdb , chainID 	= line[0][0:4], line[0][-1]
			fastaPath  		= os.path.join( fasta_dir,  'tm_%s_%s.fasta' 	% (pdb , chainID ) )
			pdbPath  		= os.path.join( pdb_dir, 	'tm_%s_%s.pdb' 		% (pdb , chainID ) )
			dict_of_pdb_objs[ pdb ].chain_List[ chainID ] = Chain( pdb, chainID, pdbPath , fastaPath, helical_fraction, tm_length)

		except KeyError:
			print error_message
			sys.exit()
	else:
			print error_message
			sys.exit()

####  pairwise comparison of all chains in each pdb, see which are redundant
# also, calc the overall helicity and size of each pdb object, store as attribute
for pdb in sorted( dict_of_pdb_objs.values() ):
	pdb.calc_size()
	pdb.calc_helicity()


	if pdb.helicity < 50 or pdb.size < 25:		# skip too short or non-helical
		beta_or_small_BADpdbs.append( pdb )
		del dict_of_pdb_objs[pdb.pdb]
		continue

	chain_group_log_txt += '%s residues = %d helicity = %.1f\n' % (pdb, int( pdb.size ), pdb.helicity)
	#print '%s residues = %d helicity = %f\n' % (pdb, int( pdb.size ), pdb.helicity)

	# parse fasta files
	pdb.parse_fasta()

	# now do each unique pairwise comparison of chains
	# For each segment, store the highest seqID between chains in matrix
	ungrouped_chains= pdb.chain_List.keys()

	for chain_1, chain_2 in combinations( sorted( pdb.chain_List.items()), 2 ):

		# skip evaluating an element that is assigned to group and is not the representative
		if chain_1[0] not in ungrouped_chains or chain_2[0] not in ungrouped_chains:
			continue

		closest_pair = {}

		for seq_1, seq_2 in product( sorted( chain_1[1].segmentDict.items()) , sorted( chain_2[1].segmentDict.items() ) ):
			
			## use blosom for global alignment... then score by % sequence identity of shorter one 
			alignment = pairwise2.align.globalds( seq_1[1] , seq_2[1], subsMat, -10, -2 , one_alignment_only=True)[0]
			seqID = 0.0
			reference_length = np.min( [len(alignment[0]), len(alignment[1])] )
			reference_length_scale 	= 1.0/reference_length
			for a, b in zip(alignment[0], alignment[1]):
				if a==b:
					seqID += reference_length_scale

			# record data just to have... 
			average_pairwise_seqID.append( seqID )

			# record the 
			try:
				if closest_pair[ seq_1[0] ][0] < seqID:
					closest_pair[ seq_1[0] ] = ( seqID, seq_2[0], reference_length )

			except KeyError:
				closest_pair[ seq_1[0] ] = ( seqID, seq_2[0], reference_length )
		

		overall_seqID = 0
		for k, v in closest_pair.items():
			overall_seqID += v[0] * v[-1]/chain_1[1].ch_size


		# combine chains to group if necessary and prune from existing list
		# representative is the longer chain, or not at all. 
		if overall_seqID > 0.92:

			# if existing chain group is taken over by new representative, 
			# make sure to transfer members to the grouping under new repr
			if chain_2[1].ch_size > chain_1[1].ch_size:
				ungrouped_chains.remove( chain_1[0] )
				pdb.chainGroups[ chain_2[0] ].append( chain_1[0] )
				try:	
					pdb.chainGroups[ chain_2[0] ].extend( pdb.chainGroups[ chain_1[0] ] )
					del pdb.chainGroups[ chain_1[0] ]
				except KeyError:
					pass

			else:
				ungrouped_chains.remove( chain_2[0] )
				pdb.chainGroups[ chain_1[0] ].append( chain_2[0] )

	# add ungrouped chains as their own groups 
	for ch in ungrouped_chains:
		pdb.chainGroups[ch].append(ch)

	# summary for PDB
	for k, v in pdb.chainGroups.items():
		#print '%s: %s' %  ( k, ' '.join( sorted( v)  ) )
		chain_group_log_txt += '%s: %s\n' %  ( k, ' '.join( sorted( v)  ) )
	chain_group_log_txt += '\n'
#	print

average_pairwise_seqID 	= np.array(average_pairwise_seqID)
len_seqID_data 			= len( average_pairwise_seqID )
med_seqID 				= round( np.median( average_pairwise_seqID ) , 3 )
mean_seqID				= round( np.mean( average_pairwise_seqID ) , 3 )
std_seqID  				= round( np.std( average_pairwise_seqID ) , 3 )


chain_group_log_txt += 'median seq_id between pairs of helices n=%d, median = %.3f  mean = %.3f +/- %.3f (std)' % (len_seqID_data, med_seqID, mean_seqID, std_seqID)

chain_group_file = open( 'good_helical_redundant_chainGroups.txt' , 'w' )
chain_group_file.write( chain_group_log_txt ) 
chain_group_file.close()

# save pickle of dict of all pdb groups with all their chain group data
cPkl__pdbData_outFile = open( 'pdb_chainGroup_data.cPickle', 'wb' )
pkl.dump( dict_of_pdb_objs, cPkl__pdbData_outFile )



