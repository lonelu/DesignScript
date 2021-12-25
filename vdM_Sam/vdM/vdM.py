import pickle

from metalprot.apps.hull import write2pymol 

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

#------------------------------
# Draw vdM score from Nick's database

import pickle
import os
import pandas as pd
from prody import writePDB, AtomGroup
import matplotlib.pyplot as plt
import numpy as np
pd.set_option("display.max_columns", None)


path_to_database = '/mnt/e/DesignData/vdM_database_20210617/'


#My trying to extract vdm score from scoring database.

path_to_scoring_list = [path_to_database, 'scoring/', 'conh2/', 'all_contacts/', 'dataframes/', 'ASP.pkl']

path_to_scoring = ''.join(path_to_scoring_list)

scoring_db = pd.read_pickle(path_to_scoring)

scores = scoring_db['C_score_bb_ind']

len(scores)

scores_unique = np.unique(scores)

with open('/mnt/e/conh2_ASP_C_score_bb_ind_score.txt', 'w') as f:
    for s in scores_unique:
        f.write(str(s)  + '\n')

# Sam's method. It extract the same infomation. 

path_to_info_list = [path_to_database, 'info/', 'conh2/', 'all_contacts/', 'dataframes/', 'ASP.pkl']
file = ''.join(path_to_info_list)
df = pd.read_pickle("%s" % file) #read pickled dataframe (*.pkl) in pandas
df1 = df.drop_duplicates(subset=['cluster_number']) # drops all redundant instances of the same value of "cluster_number" so only count each cluster_size once
Mean = df1['cluster_size'].mean() # determines the mean cluster_size of non-redundant clusters
df.insert(0, 'C', np.log((df['cluster_size'])/Mean))
df2=df
C = df.drop_duplicates(subset=['cluster_number'])[['cluster_number', 'C']]
C.to_csv('test.csv')

