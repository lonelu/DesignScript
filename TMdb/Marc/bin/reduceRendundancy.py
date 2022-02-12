import os, sys, numpy as np, cPickle as pkl, time
from collections import defaultdict
from itertools import combinations, product
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist


from tm_pdb_utilities import *

subsMat = matlist.blosum62

## NOTE:  findChainGroups.py must be put into the python path

# Look the pairwise sequence identity of PDBs TM regions, 
# by the pairwise seqID of each TM helix
# (>90% of any TM segment, or >50% overall average pairwise idenitity)

# input (1) cPickle file of all all the valid helical PDB objects created in chain groups function

# output(1) list the non-redundant representative PDBs 
# output(2)	tm pdb for each non-redundant chain, including all residue segments < 15 A away. (valid for MASTER searching) 
# output(3) tm pdb for each entire pdb, including all chains
# output(4) log file of each pdb by pdb comparison, w/ per-chain breakdown of seqID
# 			log will also note groups, representatives, and history tracking comparisons

# example command line
# python bin/reduceRendundancy.py pdb_chainGroup_data.cPickle
# WRITTEN FOR PYTHON 2.7

cPkl_file_path 	 	 = sys.argv[1]

levelSeqID 			= float( sys.argv[2] )

str_ID = '%0.f' % levelSeqID * 100

perChain_seqID_no100 = []

# load data, a directory of Pdb objects containing chain objects
dict_of_pdb_objs = pkl.load( open( cPkl_file_path , 'rb') )

ungrouped_pdbs	 = dict_of_pdb_objs.keys()
redundant_pdbs   = []

# can do it by chain as well as by whole PDB... 
ungrouped_chains = dict_of_pdb_objs.keys()

log_file_byPDB_byChain = ''
log_file_showHierarchy = ''


index = 0

for pdb_1, pdb_2 in combinations( sorted( dict_of_pdb_objs.values(), key=lambda x : x.pdb ), 2 ):


	pdb_1_chains = sorted( [ pdb_1.chain_List[x] for x in pdb_1.chainGroups.keys() ] )
	pdb_2_chains = sorted( [ pdb_2.chain_List[x] for x in pdb_2.chainGroups.keys() ] )

	redundant_flg = 0
	redundant_ch  = []

	

	if pdb_1.pdb in redundant_pdbs:
		continue
	if pdb_2.pdb in redundant_pdbs:
		continue

#	print '%s %d %.1f %s' % (pdb_1, int( pdb_1.size ), pdb_1.helicity, ' '.join( pdb_1.chainGroups.keys() ) )
#	print '%s %d %.1f %s' % (pdb_2, int( pdb_2.size ), pdb_2.helicity, ' '.join( pdb_2.chainGroups.keys() ) )
	log_file_byPDB_byChain += '%s %d %.1f %s\n' % (pdb_1, int( pdb_1.size ), pdb_1.helicity, ' '.join( pdb_1.chainGroups.keys() ) )
	log_file_byPDB_byChain += '%s %d %.1f %s\n' % (pdb_2, int( pdb_2.size ), pdb_2.helicity, ' '.join( pdb_2.chainGroups.keys() ) )

	timeRef = time.time()

	for chain_1 , chain_2 in product( pdb_1_chains, pdb_2_chains ):

		closest_pair = {}
		

		for seq_1, seq_2 in product( sorted( chain_1.segmentDict.items()) , sorted( chain_2.segmentDict.items() ) ):

			## use blosom for global alignment... then score by % sequence identity of shorter one 
			alignment = pairwise2.align.globalds( seq_1[1] , seq_2[1], subsMat, -10, -2 , one_alignment_only=True)[0]
			seqID = 0.0
			reference_length = np.min( [ len( seq_1[1] ), len( seq_2[1] ) ] )
			reference_length_scale 	= 1.0/reference_length
			for a, b in zip(alignment[0], alignment[1]):
				if a==b:
					seqID += reference_length_scale

			# record the best matching pair of helices for each helix
			try:
				if closest_pair[ seq_1[0] ][0] < seqID:
					closest_pair[ seq_1[0] ] = ( seqID, seq_2[0], reference_length )

			except KeyError:
				closest_pair[ seq_1[0] ] = ( seqID, seq_2[0], reference_length )

			if seqID > .90:
#				print '  ', chain_1 , chain_2, '\n\tMATCHING HELICES %.3f\n\t  %s\n\t  %s\n' % ( seqID, alignment[0], alignment[1] ) 
				redundant_flg = 1
				if (chain_1 , chain_2)  not in redundant_ch:
					redundant_ch.append( (chain_1 , chain_2) )

			else: 
				perChain_seqID_no100.append( seqID )
				
		

#			print seq_1, seq_2, seqID
		# derive the overall sequence identity
		overall_seqID = 0
		for k, v in closest_pair.items():
			overall_seqID += v[0] * v[-1]/chain_1.ch_size

		overall_seqID = round( 100*overall_seqID, 1 )
		perChain_seqID_no100.append( overall_seqID )

		if overall_seqID > levelSeqID * 100:
			redundant_flg = 1
			if (chain_1 , chain_2)  not in redundant_ch:
				redundant_ch.append( (chain_1 , chain_2) )

			log_file_byPDB_byChain += '\t%s %s %.1f ********\n' % ( chain_1, chain_2, overall_seqID )
		else:
			log_file_byPDB_byChain += '\t%s %s %.1f\n' % ( chain_1, chain_2, overall_seqID )


#		print '\t%s %s %.1f' % ( chain_1, chain_2, overall_seqID )

	log_file_byPDB_byChain += '\n'


	# if redundant species found, select one... lower resolution... 
	# also find any independent chains not related
	if redundant_flg:
		#print redundant_pdbs

		if pdb_1.resolution < pdb_2.resolution:
			#ungrouped_pdbs.remove( pdb_1.pdb )
			redundant_pdbs.append( pdb_2.pdb )
			log_file_showHierarchy += '%s > %s | %s\n' % ( pdb_1, pdb_2, ' '.join( [ '%s~%s' % n for n in redundant_ch ]) )
			print '%s > %s | %s \n' % ( pdb_1, pdb_2, ' '.join( [ '%s~%s' % n for n in redundant_ch ]) )
	
		else:
			#ungrouped_pdbs.remove( pdb_2.pdb )
			redundant_pdbs.append( pdb_1.pdb )
			log_file_showHierarchy += '%s > %s | %s\n' % ( pdb_2, pdb_1, ' '.join( [ '%s~%s' % n for n in redundant_ch ]) )
			print '%s > %s | %s \n' % ( pdb_2, pdb_1, ' '.join( [ '%s~%s' % n for n in redundant_ch ]) )
	
	if index in np.arange(1000, 50000, 1000):
		print index,  '\n\n'


	index +=1

log_file_Hierarchy = open( '%s_redundancyRemovalLog.txt' % str_ID  , 'w' )
log_file_Hierarchy.write( log_file_showHierarchy  )
log_file_Hierarchy.close()


log_file_byPDB = open( '%s_pdb_seq_comparisonLog.txt' % str_ID , 'w' )
log_file_byPDB.write( log_file_byPDB_byChain  )
log_file_byPDB.close()

non_redundant_txt = ''
non_redundant = []
for pdb in dict_of_pdb_objs.values():
	if pdb.pdb not in redundant_pdbs:
		non_redundant.append( pdb  )
		non_redundant_txt += '%s\n' % pdb

log_finalPdbs = open( '%s_nonRedundant_finalPdbs.txt' % str_ID , 'w' )
log_finalPdbs.write(non_redundant_txt)

print len( non_redundant ), 'non-redundant in set'

mean_seqID 		= round( np.mean( perChain_seqID_no100 ), 4 ) 
median_seqID 	= round( np.median( perChain_seqID_no100 ), 4 )
stdev_seqID 	= round( np.std( perChain_seqID_no100 ), 4 )

print 'perTMsegment seqID (if < .90), median = %.4f mean = %.4f std = %.4f' % (median_seqID, median_seqID, stdev_seqID  )






















