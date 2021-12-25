## "vdms_toPDB.py" converts the representative iFG - protein clusters into PDB (*.pdb) files that can be opened in PyMOL.
## The script should be run as: "python vdms_toPDB.py XXX.pkl" where "XXX.pkl" is the name of the pickle file you want to iterate over.
## The output filename corresponds to the following notation: "(amino acid)_(cluster size)_(log-odds score)_(vdM ID)", where the log-odds score ("C") is defined
## in "A defined structural unit enables de novo design of small-molecule–binding proteins", N.F. POLIZZI & W.F. DEGRADO, Science, 2020, 1227-1233 (see Figure 1E and SI).
import pandas as pd
import numpy as np
import sys
import math

file=sys.argv[1]
file_AA=file[0:3]

df = pd.read_pickle("%s" % file) #read pickled dataframe (*.pkl) in pandas

df1 = df.drop_duplicates(subset=['cluster_number']) # drops all redundant instances of the same value of "cluster_number" so only count each cluster_size once
Mean = df1['cluster_size'].mean() # determines the mean cluster_size of non-redundant clusters

df.insert(0, 'C', np.log((df['cluster_size'])/Mean))
df2=df

C = df2['C']

df3 = df2.astype({'iFG_count':'int'}) #converts the object to an integer for rapid grabbing of data
unique_clusters = df3["iFG_count"].unique() #saves the index of all idenfitiers based on the unique cluster_number

for x in unique_clusters:
    iFG_ID = x
    df4 = df3[df3['iFG_count'] == x]
    unique_vdMs = df4["vdM_count"].unique()
    for y in unique_vdMs:
        vdM_ID = y #saves "iFG_count" as identifier for a particular vdM for self-reference
        df5 = df4[df4['vdM_count'] == y]
        unique_AA = df5["query_name"].unique()
        for z in unique_AA:
            vdM_AA = z
            if vdM_AA == 'histidineND1':
                vdM_aa = 'HIS'
            elif vdM_AA == 'histidineNE2':
                vdM_aa = 'HIS'
            elif vdM_AA == 'tryptophan':
                vdM_aa = 'TRP'
            df6 = df5[df5['query_name'] == z]
            df6.insert(0, 'atom#', range(1,1+len(df6.index)))
            df6.insert(0, 'atom', 'ATOM')
            df6['atom_element'] = df6['name'].astype(str).str[0]  #writes atom from atom type using 1st character of object string (.str[0]) as a new column ("atom_element")
            df6.insert(2,'occupancy', 1.00) #inserts values as "float"; use '1.00' to get object datatype
            df6.insert(3,'T-factor', 0.00)
            df7 = df6[['atom','atom#','name','resname_vdm','chain','resnum','c_x','c_y','c_z','occupancy','T-factor','segment','atom_element']]
            df7_rename=df7.replace({'name':'ND1'}, 'NE2')
            df7_rename2=df7_rename.replace({'name':'NE1'}, 'NE2')
            if file_AA == 'ASP':
                df8=df7_rename2.loc[df7_rename2['name'].isin(['CG','OD1','OD2','NE2'])]
                index_names = df8[(df8['chain'] == 'Y') & (df['name'] == 'CG')].index
                df8.drop(index_names, inplace=True)
                df9=df8[['c_x','c_y','c_z']] # grabs (x,y,z) coords of each atom and creates new dataframe
                df9.index=range(4) #re-index the rows to 0 to 3 corresponding to CD, OE1, OE2, ND1
                OC1 = df9.iloc[1]-df9.iloc[0]
                unit_OC1 = OC1/np.linalg.norm(OC1)
                OC2 = df9.iloc[2]-df9.iloc[0]
                unit_OC2 = OC2/np.linalg.norm(OC2)
                NC = df9.iloc[3]-df9.iloc[0]
                unit_NC = NC/np.linalg.norm(NC)
                dot_prod_OC1NC = np.dot(unit_OC1, unit_NC)
                angle_OD1_CG_NE2 = np.arccos(dot_prod_OC1NC)
                ang_OD1_CG_NE2=math.degrees(angle_OD1_CG_NE2)
                dot_prod_OC2NC = np.dot(unit_OC2, unit_NC)
                angle_OD2_CG_NE2 = np.arccos(dot_prod_OC2NC)
                ang_OD2_CG_NE2=math.degrees(angle_OD2_CG_NE2)
                distO1N = np.sqrt(np.sum(np.power((df9.iloc[3]-df9.iloc[1]),2))) #distance between N and OE
                distO2N = np.sqrt(np.sum(np.power((df9.iloc[3]-df9.iloc[2]),2)))
                distCN = np.sqrt(np.sum(np.power((df9.iloc[3]-df9.iloc[0]),2)))
                if (distO1N < 3.0) and (distO2N < 3.0) and (distCN < 3.4) and (ang_OD2_CG_NE2 < 75) and (ang_OD1_CG_NE2 < 75):
                    fname = iFG_ID
                    fname2 = vdM_ID
                    clust_size = df5['cluster_size'].iloc[0]
                    score = df6['C'].iloc[0]
                    np.savetxt("./output_headon/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), df7, fmt=('%s','%6s','%3s','%4s','%1s','%3s','%11.3f','%7.3f','%7.3f','%5.2f','%5.2f','%6s','%4s')) # keep the "fmt" string to maintain PDB format correctly.
                    filenames = ["./helix_query.pdb", "./output/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2)]
                    with open("./helix_headon/%s_%s_%s_C=%.3f_%s_%s_helix.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), 'w') as outfile:
                        for names in filenames:
                            with open(names) as infile:
                                outfile.write(infile.read())
                            outfile.write("\n")
                if (distO1N < 3.0) and (distO2N < 5.0) and (distO2N > 4.2) and (ang_OD2_CG_NE2 > 140) and (ang_OD1_CG_NE2 > 25) and (ang_OD1_CG_NE2 < 45):
                    fname = iFG_ID
                    fname2 = vdM_ID
                    clust_size = df5['cluster_size'].iloc[0]
                    score = df6['C'].iloc[0]
                    np.savetxt("./output_sideon/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), df7, fmt=('%s','%6s','%3s','%4s','%1s','%3s','%11.3f','%7.3f','%7.3f','%5.2f','%5.2f','%6s','%4s')) # keep the "fmt" string to maintain PDB format correctly.
                    filenames = ["./helix_query.pdb", "./output/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2)]
                    with open("./helix_sideon/%s_%s_%s_C=%.3f_%s_%s_helix.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), 'w') as outfile:
                        for names in filenames:
                            with open(names) as infile:
                                outfile.write(infile.read())
                            outfile.write("\n")
                if (distO2N < 3.0) and (distO1N < 5.0) and (distO1N > 4.2) and (ang_OD1_CG_NE2 > 140) and (ang_OD2_CG_NE2 > 25) and (ang_OD2_CG_NE2 < 45):
                    fname = iFG_ID
                    fname2 = vdM_ID
                    clust_size = df5['cluster_size'].iloc[0]
                    score = df6['C'].iloc[0]
                    np.savetxt("./output_sideon/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), df7, fmt=('%s','%6s','%3s','%4s','%1s','%3s','%11.3f','%7.3f','%7.3f','%5.2f','%5.2f','%6s','%4s')) # keep the "fmt" string to maintain PDB format correctly.
                    filenames = ["./helix_query.pdb", "./output/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2)]
                    with open("./helix_sideon/%s_%s_%s_C=%.3f_%s_%s_helix.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), 'w') as outfile:
                        for names in filenames:
                            with open(names) as infile:
                                outfile.write(infile.read())
                            outfile.write("\n")
            elif file_AA == 'GLU':
                df8=df7_rename2.loc[df7_rename2['name'].isin(['CD','OE1','OE2','NE2'])]
                df9=df8[['c_x','c_y','c_z']] # grabs (x,y,z) coords of each atom and creates new dataframe
                df9.index=range(4) #re-index the rows to 0 to 3 corresponding to CD, OE1, OE2, ND1
                OC1 = df9.iloc[1]-df9.iloc[0]
                unit_OC1 = OC1/np.linalg.norm(OC1)
                OC2 = df9.iloc[2]-df9.iloc[0]
                unit_OC2 = OC2/np.linalg.norm(OC2)
                NC = df9.iloc[3]-df9.iloc[0]
                unit_NC = NC/np.linalg.norm(NC)
                dot_prod_OC1NC = np.dot(unit_OC1, unit_NC)
                angle_OE1_CD_NE2 = np.arccos(dot_prod_OC1NC)
                ang_OE1_CD_NE2=math.degrees(angle_OE1_CD_NE2)
                dot_prod_OC2NC = np.dot(unit_OC2, unit_NC)
                angle_OE2_CD_NE2 = np.arccos(dot_prod_OC2NC)
                ang_OE2_CD_NE2=math.degrees(angle_OE2_CD_NE2)
                distO1N = np.sqrt(np.sum(np.power((df9.iloc[3]-df9.iloc[1]),2))) #distance between N and OE
                distO2N = np.sqrt(np.sum(np.power((df9.iloc[3]-df9.iloc[2]),2)))
                distCN = np.sqrt(np.sum(np.power((df9.iloc[3]-df9.iloc[0]),2)))
                if (distO1N < 3.0) and (distO2N < 3.0) and (distCN < 3.4) and (ang_OE2_CD_NE2 < 75) and (ang_OE1_CD_NE2 < 75):
                    fname = iFG_ID
                    fname2 = vdM_ID
                    clust_size = df5['cluster_size'].iloc[0]
                    score = df6['C'].iloc[0]
                    np.savetxt("./output_headon/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), df7, fmt=('%s','%6s','%3s','%4s','%1s','%3s','%11.3f','%7.3f','%7.3f','%5.2f','%5.2f','%6s','%4s')) # keep the "fmt" string to maintain PDB format correctly.
                    filenames = ["./helix_query.pdb", "./output/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2)]
                    with open("./helix_headon/%s_%s_%s_C=%.3f_%s_%s_helix.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), 'w') as outfile:
                        for names in filenames:
                            with open(names) as infile:
                                outfile.write(infile.read())
                            outfile.write("\n")
                if (distO1N < 3.0) and (distO2N < 5.0) and (distO2N > 4.2) and (ang_OE2_CD_NE2 > 140) and (ang_OE1_CD_NE2 > 25) and (ang_OE1_CD_NE2 < 45):
                    fname = iFG_ID
                    fname2 = vdM_ID
                    clust_size = df5['cluster_size'].iloc[0]
                    score = df6['C'].iloc[0]
                    np.savetxt("./output_sideon/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), df7, fmt=('%s','%6s','%3s','%4s','%1s','%3s','%11.3f','%7.3f','%7.3f','%5.2f','%5.2f','%6s','%4s')) # keep the "fmt" string to maintain PDB format correctly.
                    filenames = ["./helix_query.pdb", "./output/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2)]
                    with open("./helix_sideon/%s_%s_%s_C=%.3f_%s_%s_helix.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), 'w') as outfile:
                        for names in filenames:
                            with open(names) as infile:
                                outfile.write(infile.read())
                            outfile.write("\n")
                if (distO2N < 3.0) and (distO1N < 5.0) and (distO1N > 4.2) and (ang_OE1_CD_NE2 > 140) and (ang_OE2_CD_NE2 > 25) and (ang_OE2_CD_NE2 < 45) :
                    fname = iFG_ID
                    fname2 = vdM_ID
                    clust_size = df5['cluster_size'].iloc[0]
                    score = df6['C'].iloc[0]
                    np.savetxt("./output_sideon/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), df7, fmt=('%s','%6s','%3s','%4s','%1s','%3s','%11.3f','%7.3f','%7.3f','%5.2f','%5.2f','%6s','%4s')) # keep the "fmt" string to maintain PDB format correctly.
                    filenames = ["./helix_query.pdb", "./output/%s_%s_%s_C=%.3f_%s_%s.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2)]
                    with open("./helix_sideon/%s_%s_%s_C=%.3f_%s_%s_helix.pdb" % (file_AA, vdM_aa, clust_size, score, fname, fname2), 'w') as outfile:
                        for names in filenames:
                            with open(names) as infile:
                                outfile.write(infile.read())
                            outfile.write("\n")