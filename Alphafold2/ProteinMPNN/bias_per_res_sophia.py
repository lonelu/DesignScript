pdbname = 'test'

# MPNN reads AA by mapping its index in this mpnn_alphabet string
mpnn_alphabet = "ACDEFGHIKLMNPQRSTVWY"
mpnn_alphabet_dict = {}
for ind, one_letter_AA in enumerate(mpnn_alphabet):
    mpnn_alphabet_dict[one_letter_AA] = ind

bias_by_res_pdb = {}
bias_by_res_pdb[pdbname] = {}
bias_by_res_dict = bias_by_res_pdb[pdbname]

for chain in chains:
    bias_per_residue = np.zeros([int(chain_length), 21])
    for bias_res, bias_aa in zip(bias_resnums, bias_aas):
        # because you need the index, you need the residue-1
        bias_per_residue[bias_res - 1, mpnn_alphabet_dict[bias_aa]] = 4

    bias_by_res_dict[chain] = bias_per_residue.tolist()


with open(F'{output_dir}/bias_per_res.jsonl', 'w') as f:
    json.dump(bias_by_res_pdb, f)