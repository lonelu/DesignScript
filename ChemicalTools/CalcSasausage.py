import CalcSasa
from pprint import pprint

pdb_path = '/Users/student/g-p-scaffold.pdb'
residues_list = [['X', 6], ['X', 8], ['X', 9]]

sasa_dict = CalcSasa.calc_sasa(pdb_path, 3, residues_list)
pprint(sasa_dict)
