from Bio import SeqIO

workdir = '/mnt/e/GitHub_Design/DesignScript/data/'
file = 'uniprot-reviewed-human-canonical-20200228.fasta'

count = 0
id_seq_dict = {}
with open(workdir + file) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if count < 10:
            print(record.id)
            print(record.seq)
        count += 1
        id = record.id.split('|')[1]
        id_seq_dict[id] = record.seq


with open(workdir + 'human.tsv', 'w') as f:
    for s in id_seq_dict.keys():
        f.write(s + '\t' + str(id_seq_dict[s]) + '\n')

evs = []
with open(workdir + 'evidence.txt', 'r') as f:
    count = 0
    for ev in f.readlines():
        # print(ev)
        # print(ev[:-1])
        # if count > 10:
        #     break
        if count ==0:         
            evs.append(ev[:-1] + '\t' + 'prot_seq\n')
            count +=1
            continue
        
        ev_split = ev[:-1].split('\t')
        #print(ev_split)
        if len(ev_split) > 8:
            id = ev_split[8]

            if id in id_seq_dict.keys():
                seq = str(id_seq_dict[id])
            else:
                seq = ''
            evs.append(ev[:-1] + '\t' + seq + '\n')
        count +=1

with open(workdir + 'evidence_seq.txt', 'w') as f:
    for ev in evs:
        f.write(ev)