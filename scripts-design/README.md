**Download the Rosetta Match output PDB files to your local directory (optional).**

Execute this command in a terminal on your computer.
Make staples across chain A and chain B
> python sftp_best_decoys.py /home/NetID/protein_stapling/CPG2_pAaF/CPG2-AB_pAaF-product -u NetID -p password -s 0

or make staples within chain A
> python sftp_best_decoys.py /home/NetID/protein_stapling/CPG2_pAaF/CPG2-A_pAaF-product -u NetID -p password -s 0

**Run FastDesign**

Make staples across chain A and chain B
> ../scripts/fast_design_matches.sh -pos ../CPG2/CPG2-AB -linker ../pAaF/pAaF-product -cst_suffix duplicated -nbh 8.0 -n 50 -mem 2000

or make staples within chain A
> ../scripts/fast_design_matches.sh -pos ../CPG2/CPG2-A -linker ../pAaF/pAaF-product -cst_suffix design -nbh 8.0 -n 50 -mem 2000

-nbh 8.0 is optional if your protein is small.

Wait for the calculations to be finished.

**Generate a .xls file containing the scores of the designed variant and the variants in which one designed point mutation is reverted back to the wild-type.**

Make staples across chain A and chain B
> python ../scripts/generate_scores_table.py CPG2-AB_pAaF-product -s 1 -dup

or make staples within chain A
> python ../scripts/generate_scores_table.py CPG2-A_pAaF-product -s 1

**Download the Rosetta FastDesign output PDB files to your local directory.**

Execute this command in a terminal on your computer.
Make staples across chain A and chain B
> python sftp_best_decoys.py /home/NetID/protein_stapling/CPG2_pAaF/CPG2-AB_pAaF-product -u NetID -p password -s 1

or make staples within chain A
> python sftp_best_decoys.py /home/NetID/protein_stapling/CPG2_pAaF/CPG2-A_pAaF-product -u NetID -p password -s 1

**Take a look of the Rosetta-designed models, manually revert those irrelevant mutations.**

1) Take a look at the overall scores and enzdes constraint scores in the .xls file. Rule out those designs with high overall scores or enzdes constraint scores.
2) For those deisgns with low scores, open the output PDB files in PyMOL one by one. Compare them with the wild-type protein sequence. Find those unnecessary mutations far away from the UAA staple.
3) Remove those unnecessary mutations in the manual.sh shell script. Run the modified shell script one by one.
> cd CPG2-A_pCaaF-product/X*Z*
> vi manual.sh
> bash manual.sh

Wait for the calculations to be finished.

**Add manually designed results to the existing .xls file.**

Make staples across chain A and chain B
> python ../scripts/generate_scores_table.py CPG2-AB_pAaF-product -s 2

or make staples within chain A
> python ../scripts/generate_scores_table.py CPG2-A_pAaF-product -s 2

**Download the manually designed output PDB files to your local directory.**

Execute this command in a terminal on your computer.
Make staples across chain A and chain B
> python sftp_best_decoys.py /home/NetID/protein_stapling/CPG2_pAaF/CPG2-AB_pAaF-product -u NetID -p password -s 2

or make staples within chain A
> python sftp_best_decoys.py /home/NetID/protein_stapling/CPG2_pAaF/CPG2-A_pAaF-product -u NetID -p password -s 2
