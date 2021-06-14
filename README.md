#  Rosetta Match and FastDesign pipeline
## folder structure ##
    /scripts-match
        generate_position_files.py
        match.sh
        general_match.flags
        generate_fast_design_input.sh
        generate_fast_design_input.py
    /scripts-design
        find_match_intersection.py
        fast_design_matches.sh
        fast_design.py
        revert_designed_residues.sh
        get_the_best_design.py
        generate_scores_table.py

    /proteins
        /${protein}
            ${protein}.pdb
            ${protein}.symm // Optional, for symmetric proteins.
            ${protein}-A // Position files for matching within chain A. Created by generate_position_files.py
            ${protein}-AB // Position files for matching across chain A and chain B. Created by generate_position_files.py

    /ligands
        /${ligand}
            /${ligand}-intermediate
                ...
            /${ligand}-TS
                ...
            /${ligand}-product
                subs.flags
                ${ligand-3-letter-name}_ligand.params
                ${ligand-3-letter-name}.rotlib.pdb
                ${matching-residue1-3-letter-name}.params
                ${matching-residue2-3-letter-name}.params
                ${ligand}-product.cst
                ${ligand}-product_duplicated.cst // Optional, for symmetric proteins.

    /designs // this folder name can be customized
        /${protein-scaffold}_${ligand} // manually created
            /${protein-scaffold}-A_${ligand}-intermediate // Created by match.sh.
                ...
            /${protein-scaffold}-A_${ligand}-TS // Created by match.sh.
                ...
            /${protein-scaffold}-A_${ligand}-product // Created by match.sh.
                /${protein-scaffold}-A_${n} // n=1,2,3... Created by match.sh Contains Match output CloudPDB files.
                ...

                /X${position1}Z${position2} // Created by generate_fast_design_input.sh. Contains FastDesign input files of each match.
                    X${position1}Z${position2}.pdb // Created by generate_fast_design_input.sh.

                    ${ligand-3-letter-name}.rotlib.pdb // Created by generate_fast_design_input.sh.

                    /design // Created by fast_design_matches.sh. Contains FastDesign output files for the match.

                    ${protein-scaffold}-A_${point_mutation1}_${point_mutation2}_${point_mutation3}_${point_mutation4}.pdb
                        // Created by revert_designed_residues.sh when calling get_the_best_design.py. The best decoy of all FastDesign output decoys.

                    /revert_${point_mutation3} // Revert each designed point mutations back to the wild type residue. Created by revert_designed_residues.sh.
                    /revert_${point_mutation4}
                    ...
                /X${position5}Z${position6}
                    ...
                ...
            /${protein-scaffold}-AB_${ligand}-intermediate // created by match.sh
                ...
            /${protein-scaffold}-AB_${ligand}-TS // created by match.sh
                ...
            /${protein-scaffold}-AB_${ligand}-product // created by match.sh
                ...
