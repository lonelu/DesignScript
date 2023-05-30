Install Phenix and Coot (WinCoot in windows.)

1. run xtriage for data quanlity using aimless.mtz.
    a. get Copy Number from 'Solvent content and Matthews coefficient'. We want the Solvent conctent around 0.5.
    b. The I/sigI>15 > 1 tells the potential resolution. I/sigI>2 should be > 80.

2. run phaser for molecula replacement.
    The Copies to search for need to be changed to the Copy number obtained from xtriage.
    a. the Top TFZ should >8.

3. run autobuild. (There is no need to run autobuild if MR)
    Uisng the overall_best_placed.pdb and overall_best_denmod_map_coeffs.mtz

4. run refine a few times. 
    Each time change the pdb file to the new generated one.
    No need to select add H and update water. (Remove H and H2O if added. Do this in the last step.)

5. run coot refinement.
    a. Draw cell-sym on Unit Cell.
    b. Using GoToAtomRefine to fix from 1st to last.
    c. Measure->PointDistOn to add H2O.
    d. Calc -> OtherModeling -> FindH2O -> AddNew
    e. PlaceAtomAs -> SO4 -> Add.

    Remove red bar based on the 3 plots. 

6. run refine a few times.
    In the refine result, there is indications for clashes and impect things which can be refered coot.

7. run coot. Add H2O and salt finally.