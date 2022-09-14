From Sophia

I don’t have a manual, but after you load your ligand, on the menu, you can go to Protein -> Structure Prepation.
In the Structure Preparation window, you can click “Protonate3D” so that it will add missing hydrogens, and calculate partial charges. (In the Protonate 3D menu, I just use all the defaults).  If you get an error in the Structure Preparation like “non-unique name”, then it means MOE put in atoms that have the same name as another atom, so you can just click that line and then click “Correct”
Sometimes it gets the formal charges or the hydrogens wrong, so after Protonate 3D, you should look at your ligand to make sure the hydrogens are placed correctly. When you click on each atom, on the bottom of the screen, it’ll tell you the atom name and atom type. the atom type will show a charge if the atom is charged. I think for your ligand, there should be no charge.
to get the conformers, on the menu, go to Compute -> Conformations -> Search.  I use all the defaults, and then I just run it. You can try it with and without QM.

oh yeah if you want you can change the forcefield (bottom left corner) to Amber14:EHT. it uses amber14 for protein, EHT for small molecule. I don’t know too much about the other forcefields.