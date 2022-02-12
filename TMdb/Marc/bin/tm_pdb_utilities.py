import os, sys, numpy as np
from collections import defaultdict

class Pdb():
		
	def __init__(self, pdb_name, exp_type, resolution, pdb_dir):
		self.pdb		= pdb_name
		self.chain_List	= {}
		self.newPDBpath = os.path.join( pdb_dir, '%s.pdb' % self.pdb )
		self.helicity	= 'NOT YET CALCULATED'
		self.size		= 'NOT YET CALCULATED'
		self.expType 	= exp_type
		self.resolution = resolution
		self.chainGroups= defaultdict(list)

	def __repr__(self):
		return self.pdb

	def calc_helicity(self):
		if isinstance(self.size, basestring) or len(self.chain_List) == 0:
			self.helicity = 0
		else:
			self.helicity = round( 100 * np.sum( [ x.ch_size * x.helix_frac * 0.01 / self.size for x in self.chain_List.values() ] ), 1)
		return self.helicity

	def calc_size(self):
		if len(self.chain_List) == 0:
			self.size = 0
		else:
			self.size = np.sum( x.ch_size for x in self.chain_List.values() )
		return int( self.size )

	# parse the fasta files for each chain. 
	# Store each helix seq in a list, hashed to each segment ID
	def parse_fasta(self):
		
		flag = 0
		for chId, obj in self.chain_List.items():
			for line in open( obj.fastaPath ):
				
				if line[:2] == '>>':		# store header
					segment_info = line[2:].rstrip()
				else:
					obj.segmentDict[ segment_info ] = line.rstrip() 

class Chain():
	def __init__(self, pdb_name, chain, pdbPath , fastaPath, helical_fraction, tm_length):
		self.pdb 		= pdb_name
		self.chain 		= chain
		self.pdbPath 	= pdbPath
		self.fastaPath	= fastaPath
		self.helix_frac = helical_fraction
		self.ch_size	= tm_length
		self.segmentDict= defaultdict(list)

	def __repr__(self):
		return '%s_%s' % ( self.pdb, self.chain )

UnNatAA={}
UnNatAA["ALA"] = 'A'; UnNatAA["CYS"] = 'C'; UnNatAA["ASP"] = 'D'; UnNatAA["GLU"] = 'E'; UnNatAA["PHE"] = 'F';
UnNatAA["GLY"] = 'G'; UnNatAA["HIS"] = 'H'; UnNatAA["ILE"] = 'I'; UnNatAA["LYS"] = 'K';
UnNatAA["LEU"] = 'L'; UnNatAA["MET"] = 'M'; UnNatAA["ASN"] = 'N'; UnNatAA["PRO"] = 'P'; UnNatAA["GLN"] = 'Q';
UnNatAA["ARG"] = 'R'; UnNatAA["SER"] = 'S'; UnNatAA["THR"] = 'T'; UnNatAA["VAL"] = 'V'; UnNatAA["TRP"] = 'W'; UnNatAA["TYR"] = 'Y';
UnNatAA['ABA'] = 'A'; UnNatAA['CSO'] = 'C'; UnNatAA['CSD'] = 'C'; UnNatAA['CME'] = 'C';
UnNatAA['OCS'] = 'C'; UnNatAA["HSD"] = 'H'; UnNatAA['KCX'] = 'K'; UnNatAA['LLP'] = 'K';
UnNatAA['MLY'] = 'K'; UnNatAA['M3L'] = 'K'; UnNatAA['MSE'] = 'M'; UnNatAA['PCA'] = 'P'; UnNatAA['HYP'] = 'P';
UnNatAA['SEP'] = 'S'; UnNatAA['TPO'] = 'T'; UnNatAA['PTR'] = 'Y'


