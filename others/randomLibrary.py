'''
For sequence mutation detection. (data from Ketaki) 
'''

import os
import pandas as pd


workdir = '/mnt/e/_temp/'
#infile = 'Random library.txt'
#outfile = 'RandamLibraryOut2.txt'
infile = 'Sort round 4.txt'
outfile = 'Sort_Round4_Out'


aa20 = 'GAVLIMFWPSTCYNQDEKRH'
std_seq = 'SEDEKLADTIKKKYAEAAAVGQEAVAVFNTMKAAFQNGDKEAVAQYLARLASLYTRHEELLNRILEDAKKAGWEKAVTLMNEFTATFQTGKSIFNAMVAAFKNGDDDSFESYLQALEKVTAKGETLADKIEKAIEEAIQKLK'


seqs = []
with open(workdir + infile, 'r') as f:
    #The header is manualy removed. Otherwise header need to be removed by skip the first line.
    for line in f.readlines():
        seq = line.split('\n')[0].split('*')[0]
        if len(seq) != len(std_seq):
            continue
        seqs.append([c for c in seq])


df = pd.DataFrame(seqs)
total = df.shape[0]

stats = []
for i in range(len(std_seq)):
    std_c = std_seq[i]
    stat = []
    for aa in aa20:
        if aa == std_c:
            #count = -total
            count = -1
        else:
            count = df.loc[df. iloc[:, i] == aa].shape[0]
        #stat.append(round(count/total, 4))
        stat.append(count)
    stats.append(stat)

df_out = pd.DataFrame(stats)
df_out = df_out.rename(index = lambda x:  std_seq[x], columns=lambda x: aa20[x])

df_out.to_csv(workdir + outfile + '_total_' + str(total) + '.txt', sep='\t')

