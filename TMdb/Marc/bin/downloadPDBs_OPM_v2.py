# Download all TM domain containing PDBs from the Orientations of Membrane Protein Database
# 03/22/20, from http://opm.phar.umich.edu/subunits.php
# http://arteni.cs.dartmouth.edu/cccp/output/e8846cb43e18.html

# python bin/downloadPDBs_OPM.py allTMsegments_OPM.txt
# WRITTEN FOR PYTHON 2.7

### MAKE SURE tm_pdb_utilities.py is in your $PYTHON_PATH
 
import sys, os, subprocess as sp, re, numpy as np, shutil
from collections import Counter, defaultdict
from prody import *

from tm_pdb_utilities import *

files_w_badSegments = ['6nbq']

# parse input file and download

logFile_text = ''

downloadDir = 'selectedOPM'
fastaDir 	= 'selectedfasta'
headerDir 	= 'selectedHeaders'
seen		= []


if not os.path.exists( downloadDir ):
	os.mkdir(downloadDir )

if not os.path.exists( fastaDir ):
	os.mkdir( fastaDir )

if not os.path.exists( headerDir ):
	os.mkdir( headerDir )

good_xray_pdbList 	= []
all_exp_data_types 	= {}
good_PDB_expType 	= {}
good_PDB_byChain	= defaultdict(list)

good_exp_types = ['ELECTRON CRYSTALLOGRAPHY', 'X-RAY DIFFRACTION', 'ELECTRON MICROSCOPY']

with open(sys.argv[1]) as fin:
	for i in fin:

		if i[0] == 'i' or len(i) < 2: continue

		line = i.split(',')
		pdb 	= line[5][4:8]
		chainID = line[6]

		# debug option
		#if pdb != '6nbq':
		#	continue

		print '\n--------------------------------\nlooking at:', pdb, chainID

		targetPath 		= os.path.join(  downloadDir, '%s.pdb' % pdb )
		fasta_address 	= 'https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=%s&compressionType=uncompressed' % pdb.upper()
		fasta_path 		= os.path.join( fastaDir, '%s.fasta' % pdb.upper() )
		header_address 	= 'https://files.rcsb.org/header/%s.pdb' % pdb.upper()
		header_path 	= os.path.join( headerDir ,'%s.txt' % pdb.upper() )
		pdb_address2 	= 'https://opm-assets.storage.googleapis.com/pdb/%s.pdb' % pdb.lower()
		pdb_address 	= '/Users/mmravic/Desktop/OPM_pdb/%s.pdb' % pdb.lower()
		pdb_path 		= os.path.join( downloadDir ,'%s.pdb' % pdb.upper() )
		good_xray_flag  = 0



		newPDBPath = os.path.join( downloadDir, 'tm_%s_%s.pdb' % (pdb.upper(), chainID.upper()) )

		# if this is not the first time running this file, skip already processed entries
		if os.path.exists(newPDBPath):
			print 'already parsed', newPDBPath, '# usable so far:',len( good_xray_pdbList ) , ' '.join(line[11:-3]) , '...'
			if pdb not in good_xray_pdbList:
				good_xray_pdbList.append(pdb)
		
		#	continue

#		print 'now looking at', newPDBPath,'\n # good PDBs so far...',  len( good_xray_pdbList ) ,'...\n'

		if not os.path.exists( header_path ):						# Download file, but only if not too large (check fasts if > 5000 residues)
			download_header_command  = ['wget', '-O', header_path, header_address ]
			sp.call( download_header_command )

		#print pdb, chainID,# good_xray_pdbList

		# check in header if structure is a good X-ray (3.3 A) or electron diffraction struct (<3.1 A)
		# store whether this has been checked or not
		if pdb.upper() in seen:
			if pdb.upper() not in good_xray_pdbList:
				continue
		else:
			seen.append( pdb.upper() )
			with open( header_path ) as headerFile:
				for line in headerFile:
					if line[:6] == 'EXPDTA':
						exp_type = line[6:].rstrip().lstrip()
		

						all_exp_data_types[pdb.upper()] = exp_type

						if exp_type not in good_exp_types:
							print '----- SKIP!!', pdb
							break 
						else:
							pass


					if line[:21] == 'REMARK   2 RESOLUTION':
						resolution = float( line.split()[-2] )
						print exp_type, resolution
						if exp_type == 'X-RAY DIFFRACTION' and resolution <= 3.3:
							
							good_xray_pdbList.append( pdb.upper() )
							good_xray_flag = 1

							all_exp_data_types[pdb.upper()] = '%s %s' % ( exp_type, str(resolution) )

						elif resolution <= 3.1:

							print 'good'
							good_xray_pdbList.append( pdb.upper() )
							good_xray_flag = 1

							all_exp_data_types[pdb.upper()] = '%s %s' % ( exp_type, str(resolution) )
						else:
							print '----- SKIP!!', pdb, '\n'

			if not good_xray_flag:
				seen.append( pdb )
				continue

		# download fasta and pdb if structure fits criteria
		if not os.path.exists( fasta_path ):
				download_fasta_command = ['wget', '-O', fasta_path, fasta_address]
				sp.call( download_fasta_command )

		if not os.path.exists( pdb_path ):
				try:
					shutil.copy( pdb_address, pdb_path )
				except IOError:
					download_pdb_command = ['wget', '-O', pdb_path, pdb_address2]
					sp.call( download_pdb_command )

		
		# use the line in the OPM TM helix definitions file to select only the TM region of each PDB 
		# write to new file in selectedPDB dir... and fasta for each segment in selectedFasta

		logFile_text += "\n\n\n%s\t%s\t%s\t%s\n"%( pdb, chainID, exp_type, resolution )
		print "\n\n\n%s\t%s\t%s\t%s\n"%( pdb, chainID, exp_type, resolution )

		# parse input definitions for regex of residue ranges, split by TM helix/segment
		initial_tm_segments = {}
		match = re.findall('\((\s*\d+)-(\s*\d+)\)', i[7:])
		if not match:
			print ('\nERROR found parsing input TM definitions file\nThe definitions should fit in this regex: \' \((\d+)-(\d+)\)\'')
			print (i,'\n')
			sys.exit()
		else:
			for n, residues in enumerate( match, 1 ):
				initial_tm_segments[n] = [ int(x) for x in list(residues) ]
				#print residues

		## find each TM domain within ProDy atom group
		# add additional 4 residues in segment at N and C termini if helical
		# NOTE: REMEMBER TO REMOVE HIGHLY KINKED HELICES FROM PAIRS
		
		try:
			in_PDB 			= parsePDB( pdb_path , chain=chainID).select('protein').copy()
		except AttributeError:

			print 'HEY BIG PROBLEM HERE... PDB was not parse... skipping but user should know about this...\n\n\nREALLY!\n\n***\n\n\n\n'
			continue

		newPDB  		= 0
		totalHelical 	= 0
		fastaFile_txt 	= ''
		segment_number  = 1

		#print in_PDB.select('ca').getResnums()
		#print in_PDB.select('ca').getResnames()

		for n, resi_range in sorted( initial_tm_segments.items(), key=lambda x: x[0] ):
			selection_str = 'resnum %d to %d' % (resi_range[0], resi_range[1] )						# define selection string, residue ranges
			N_term_resis  = np.arange( resi_range[0] - 4, resi_range[0] + 1 ) 
			C_term_resis  = np.arange( resi_range[1], resi_range[1] + 5 )

			
			segment 	= in_PDB.select( selection_str )								#Select AtomGroup objects by residue number

			### Some segments have invalid residues, like 'UNK', which are ignored by ProDy
			## e.g. 6nbq	A    ('   7', '  17');  skip these, but write note in log

			#print selection_str, resi_range, N_term_resis , C_term_resis , segment
			if segment == None:
				print "---'THIS SEGMENT IS NOT VALID - ERROR IN DEFINITIONS OR INVALID RESIDUES (i.e. UNK)\n%s\n"%selection_str
				logFile_text += "---'THIS SEGMENT IS NOT VALID - ERROR IN DEFINITIONS OR INVALID RESIDUES (i.e. UNK)\n%s\n"%selection_str
				
				# hand-check these examples - fix TM segi definitions or confirm UNK residues
				if pdb in files_w_badSegments:
					print "--- marked as okay to skip, i.e. INVALID RESIDUES (i.e. UNK)\n"
					logFile_text += "--- marked as okay to skip, i.e. INVALID RESIDUES (i.e. UNK)\n"
					continue
				else:
					sys.exit()
			else:
				pass

			# Check if terminal extension is okay to add by 3 helic dihedral angles, 
			#    turn on flag for each residue to add or not add at the end, 1st & 2nd share a flag
			# ALSO skip the residues that are negative in PDB code... since ProDy throws errors selecting these resnums
			N_term_flg = np.zeros(3)
			for n, resnum in enumerate( N_term_resis[1:-1], 1 ) :

				resi = in_PDB[chainID, resnum]


				if resi:
					try:
						phi, psi = calcPhi(resi), calcPsi(resi)
						if -120 < phi < -35 and -75 < psi < 15 and resnum > 0:
							N_term_flg[n-1] = 1
					except (TypeError, ValueError):
						pass

			C_term_flg = np.zeros(3)
			for n, resnum in enumerate( C_term_resis[1:-1], 1 ) :
				resi = in_PDB[chainID, resnum]
				if resi:
					try:
						phi, psi = calcPhi(resi), calcPsi(resi)
						if -120 < phi < -35 and -75 < psi < 15:
							C_term_flg[n-1] = 1
					except (TypeError, ValueError ):
						break

			# add extension residues into atom group... extract sequence
		 	cnt = 1
		 	newAtomGroup = 0
#		 	print N_term_flg, C_term_flg, n
			while N_term_flg[-1]:
#				print cnt, N_term_flg

				# Make new atom group starting with the extra residue(s)
				if N_term_flg[-cnt] and not newAtomGroup:
					
					# check if current residue is -3, if so include both -4 and -3; else just add it (i.e. -2 or -1)
					if cnt == 3:
						resnum = ' '.join( str(n) for n in N_term_resis[:-3] if n > 0 )
					else:
						resnum = N_term_resis[-cnt-1]
					
					newAtomGroup 	= in_PDB.select( 'resnum %s' % resnum )

				# append new residue(s) to N-termini of the existing extension
				elif N_term_flg[-cnt] and newAtomGroup:

					# check if current residue is -3, if so include both -4 and -3; else just add it (i.e. -2 or -1)
					if cnt == 3:
						resnum = ' '.join( str(n) for n in N_term_resis[:-3] if n > 0)
					else:

						resnum = N_term_resis[-cnt-1]	


					newResi  	 = in_PDB.select( 'resnum %s' % resnum )
					newAtomGroup = newResi + newAtomGroup

				elif not N_term_flg[-cnt]:
					break 

				else:
					break 

				# safety debugging note
				if cnt == 3: break

				cnt += 1

			if not newAtomGroup:
				newAtomGroup = segment
			else:
				newAtomGroup = newAtomGroup + segment


			# repeat same addition process
			cnt = 0
			while C_term_flg[0]:
#				print cnt, C_term_flg,

				# Make new atom group starting with the extra residue(s)
				if C_term_flg[cnt]:
					
					if cnt == 2:
						resnum = ' '.join( str(n) for n in C_term_resis[3:] )
					else:
						resnum = C_term_resis[cnt+1]

					newResi  	 = in_PDB.select( 'resnum %s' % resnum )
					newAtomGroup = newAtomGroup + newResi

				else:
					break 

				# safety debugging note
				if cnt ==2: break

				cnt += 1


			# record new resnum range, length, and fasta segment information
			# record percent of alpha helical residues in segment
			tmSeg_seq 			 = ''
			cnt, resRange, helical = 0.0, 0, 0

			for res in newAtomGroup.copy().iterResidues():
				tmSeg_seq += UnNatAA[res.getResname()]
				resnum	= res.getResnum()
				try:
					phi, psi = calcPhi(res), calcPsi(res)
					if -100 < phi < -35 and -75 < psi < -10:
						helical += 1
				except (TypeError, ValueError ):
						pass

				if not resRange:
					resRange = (resnum, 0)
			resRange = (resRange[0], resnum)
			totalHelical += helical
			identifier = '%s_%s ' % (pdb, chainID) + '%d-%d' % resRange
			fastaFile_txt += '>>%s\n%s\n' % ( identifier, tmSeg_seq )
			logFile_text  += '\t%s\t%s\n' % ( identifier, tmSeg_seq )
 
			# set the segment ID number for this TM segment
			newAtomGroup = newAtomGroup.copy()
			newAtomGroup.setSegnames( [str(segment_number) for x in newAtomGroup.getSegnames()] )


			if not newPDB:
				newPDB = newAtomGroup
			else:
				newPDB += newAtomGroup

			segment_number  += 1

		newPDB = newPDB.copy()
		helix_percent = 100*round( float(totalHelical) / newPDB.numResidues() , 2 )

		chain_summary = '%s_%s num_residues=%d percent_helix: %.0f' % (pdb.upper(), chainID, newPDB.numResidues(), helix_percent )
		good_PDB_byChain[pdb.upper()].append( chain_summary )
		
		title =  'tm_%s_%s %s %.1f num_residues=%d percent_helix: %.0f' % (pdb.upper(), chainID, exp_type, resolution, newPDB.numResidues(), helix_percent )
		newPDB.setTitle(title)
		#newPDBPath = os.path.join( downloadDir, 'tm_%s_%s.pdb' % (pdb.upper(), chainID.upper()) )
		newFastaP  = os.path.join( fastaDir, 'tm_%s_%s.fasta' % (pdb.upper(), chainID.upper()) )

		outFasta = open( newFastaP, 'w' )
		outFasta.write( fastaFile_txt )
		outFasta.close()
		
		writePDB( newPDBPath, newPDB )
		logFile_text  += '%s\n' % title

		#print logFile_text, '\n\n'
		

			

outLog 			= open( 'parsing_TM_pdbLOG.txt', 'w' )
outLog.write( logFile_text )

# write list of all pbds with thier experimental data type + info, 
# and another of just the acceptable ones
good_pdb_file 	= open( 'suitable_redundant_pdbSet.txt', 'w' )
good_pdb_byCh 	= open( 'suitable_redundant_pdbByChain.txt', 'w' )
all_pdb_files 	= open( 'dataTypes_pdbSet.txt', 'w' )

exp_type_txt 	= ''
good_pdb_txt	= ''
good_byChTxt 	= ''


for k, v in sorted( all_exp_data_types.items() ):


	if k in good_xray_pdbList:
		good_pdb_txt += '{:<5} {:<25}\n'.format(k,v)
		exp_type_txt += '{:<5} {:<25} Include\n'.format(k,v)
		good_byChTxt += '{:<5} {:<25}\n'.format(k,v)
		for ch_info in good_PDB_byChain[k]:
			good_byChTxt += '\t{:<50} \n'.format(ch_info)

	else:
		exp_type_txt += '{:<5} {:<25} Remove\n'.format(k,v)

all_pdb_files.write( exp_type_txt[:-1] )
good_pdb_file.write( good_pdb_txt[:-1] )
good_pdb_byCh.write( good_byChTxt[:-1] )


print len( good_xray_pdbList ), 'potentially redunant, but acceptable pdb models'
print 'TM pdbs made per chain, TM fastas made per chain... \nNeed to group redundant chains and exclude similar pdbs'
print '\nDone!\n\n\n'




### remaining things to do...
# parse







