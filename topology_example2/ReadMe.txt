#This is learned from bootcamp hosted by Sam M in 12/14/2020.

> ls
all_ala_resfile.txt backbone_ZnNDI_start_round1_57_round2_92.pdb  options_all_ala_fixbb.txt

> ~/bin/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/fixbb.static.linuxgccrelease @options_all_ala_fixbb.txt -ignore_unrecognized_res
#the -ignore_unrecognized_res ignore the ligand that is not specified for Rosetta.
#generate a new file: backbone_ZnNDI_start_round1_57_round2_92_all_ala.pdb

#copy and make changes in topology_usage.txt file, based on the restrictions, filepath, etc.
#create resfile.txt, copy pkl file, copy topology.py file.
> ipython
#inside ipython, copy code in topology_usage.txt file and run each line.

#Generate resfile_topology.txt. The file contains changable amino acids for each residue.