import pickle 

#-------------------------------
filename = '/mnt/e/GitHub_Design/vdM/Ideal_ALA/ALA_pdb.pkl' 
filename = '/mnt/e/DesignData/combs_database_/clusters/hb_only/backboneCO/20181002/PHI_PSI/ALA.pkl'

filename = '/mnt/e/DesignData/combs_database_/relative_vdms/hb_only/carboxamide/20181002/HNCA.pkl'
infile = open(filename, 'rb') 
new_dict = pickle.load(infile)
infile.close()
print(new_dict)
print(new_dict.keys())
#new_dict['name']
#-------------------------------



file = '/mnt/e/DesignData/combs_database_/clusters/hb_only/backboneCO/20181002/PHI_PSI/ALA.pkl'

file = '/mnt/e/DesignData/combs_database_/clusters/hb_only/backboneCO/20181002/PHI_PSI/ALA.pkl'

file_AA = 'ALA'