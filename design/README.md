**Run FastDesign**

Make staples across chain A and chain B
> ../scripts/fast_design_matches.sh -pos ../CPG2/CPG2-AB -linker ../pAaF/pAaF-product -symm ../CPG2/CPG2_relaxed.symm -nbh 8.0 -n 50 -mem 2000

or make staples within chain A
> ../scripts/fast_design_matches.sh -pos ../CPG2/CPG2-A -linker ../pAaF/pAaF-product -symm ../CPG2/CPG2_relaxed.symm -nbh 8.0 -n 50 -mem 2000

-nbh 8.0 is optional if your protein is samll.

**Run FastRelax with designed point mutations reverted back to the wild type one by one**

Make staples across chain A and chain B
> ../scripts/revert_designed_residues.sh -pos ../CPG2/CPG2-AB -linker ../pAaF/pAaF-product -symm ../CPG2/CPG2_relaxed.symm -nbh 8.0 -n 50 -mem 2000 -get_best true

or make staples within chain A
> ../scripts/revert_designed_residues.sh -pos ../CPG2/CPG2-A -linker ../pAaF/pAaF-product -symm ../CPG2/CPG2_relaxed.symm -nbh 8.0 -n 50 -mem 2000 -get_best true

**Generate a .xls file containing the scores of the designed variant and the variants in which one designed point mutation is reverted back to the wild type**
Make staples across chain A and chain B
> python ../scripts/generate_scores_table.py CPG2-AB_pAaF-product -params ../pAaF/pAaF-product/AAF_design.params ../pAaF/pAaF-product/CYX.params ../pAaF/pAaF-product/TYZ.params

or make staples within chain A
> python ../scripts/generate_scores_table.py CPG2-A_pAaF-product -params ../pAaF/pAaF-product/AAF_design.params ../pAaF/pAaF-product/CYX.params ../pAaF/pAaF-product/TYZ.params

**Download the output PDB files to your local directory**

Execute this command in a terminal on your computer.
Make staples across chain A and chain B
> python sftp_best_decoys.py /home/NetID/protein_stapling/CPG2_pAaF/CPG2-AB_pAaF-product -u NetID -p password

or make staples within chain A
> python sftp_best_decoys.py /home/NetID/protein_stapling/CPG2_pAaF/CPG2-A_pAaF-product -u NetID -p password
