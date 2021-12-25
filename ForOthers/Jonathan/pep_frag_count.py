
'''
# Chang the work directory, Copy the command and run it in your python enviroment.
# Change script path before you run the script.
python /mnt/e/GitHub_Design/DesignScript/Jonathan/pep_frag_count.py
'''

workdir = '/mnt/e/GitHub_Design/DesignScript/Jonathan/'  # Change here if necessary

filepath = workdir + 'hitlist_collated.txt' # Change here if necessary.

output_path = workdir + 'fragment_count3.txt' # Change here if necessary.

output_sequence = False # Change here to True if you want to print sequence in the separate raw.


### function for generating sequence fragments

def get_seq_frags(sequence, min = 3):
    fragments = []
    seq_len = len(sequence)

    for frag_len in range(min, seq_len + 1):
        for i in range(0, seq_len - frag_len + 1):
            frag = sequence[i: i + frag_len]
            fragments.append(frag)

    return fragments

'''
#Test
seqs = 'ABCDEFG'

frags = get_seq_frags(seqs)
print(frags)
# should be ['ABC', 'BCD', 'CDE', 'DEF', 'EFG', 'ABCD', 'BCDE', 'CDEF', 'DEFG', 'ABCDE', 'BCDEF', 'CDEFG', 'ABCDEF', 'BCDEFG', 'ABCDEFG']
'''

### Load the data into dictionary. Here assume all the sequence is unique. If not, the code needs to be changed.

seq_dict = {} #{ABCDEFG : 0}
with open(filepath, 'r') as f:
    count = 0
    for line in f.readlines():
        if count == 0:
            count += 1
            continue
        seq = line.split('\n')[0]
        
        seq_dict[seq] = count
        count += 1

print(len(seq_dict.keys()))

### Generate fragment for each sequence and put it into dictionary 

frag_dict = {} # {'ABC': seq_ids = [1, 2, 3, 4, 5], len(ABC) = 3, count = 5, seq_seqs = [ABCDEFG ...] } The fragment ABC is from sequence 1, 2, 3, 4, 5. The len is 3
for k in seq_dict.keys():
    v = seq_dict[k]

    frags = get_seq_frags(k)

    for fg in frags:
        if fg in frag_dict.keys():
            frag_dict[fg][0].append(v)
            frag_dict[fg][2] += 1
            frag_dict[fg][3].append(k)
        else:
            frag_dict[fg] = [[v], len(fg), 1, [k]]

### Write the fragment infomation into txt file

with open(output_path, 'w') as f:
    f.write('FragmentSequence\tSequenceId\tLength\tCount\n')
    for k in frag_dict.keys():
        v = frag_dict[k]
        f.write(k + '\t' + '||'.join([str(a) for a in v[0]]) + '\t' + str(v[1]) + '\t' + str(v[2]) + '\t' + '||'.join([a for a in v[3]]) + '\n')
        if output_sequence:
            for a in v[3]:
                f.write('\t' + a + '\n')