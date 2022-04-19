# From Sam S

gmx trjconv -s md_0_10.tpr -f md_0_10.trr -o md_0_10_center.xtc -center -pbc mol -ur compact
gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact
gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o start.pdb -dump 0 (this file gets used for a reference in the next line to load the xtc file, so ensure you grab the same atoms for the output!)
gmx trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o md_0_10_fit.xtc -fit rot+trans
LOAD PDB of start frame!!
Then use "load_traj XXXXX.xtc" to grab whole trajectory. Note must navigate to the correct directory while in Pymol.



https://manual.gromacs.org/documentation/current/onlinehelp/gmx-trjconv.html
This is how I extract the frames in a less precise format (.xtc) and use the protein for self-centering From Gromacs