# Pipeline for Rosetta Match and subsequent FastDesign
## folder structure ##
    /scripts
        generate_position_files.py
        match.sh
        general_match.flags
        generate_fast_design_input.sh
        generate_fast_design_input.py
        find_match_intersection.py
        fast_design_matches.sh
        fast_design.py
        revert_designed_residues.sh
        get_the_best_design.py
        generate_scores_table.py

    /${protein-scaffold}
        ${protein-scaffold}.pdb
        ${protein-scaffold}.symm // Optional, for symmetric proteins.
        ${protein-scaffold}-A // Position files for matching within chain A. Created by generate_position_files.py
        ${protein-scaffold}-AB // Position files for matching across chain A and chain B. Created by generate_position_files.py

    /${ligand}
        /${ligand}-intermediate
            ...
        /${ligand}-TS
            ...
        /${ligand}-product
            ${ligand-3-letter-name}.rotlib.pdb
            ${ligand-3-letter-name}.params
            ${matching-residue1-3-letter-name}.params
            ${matching-residue2-3-letter-name}.params
            subs.flags
            ${ligand}-product.cst
            ${ligand-3-letter-name}_design.params
            ${ligand}-product_design.cst // Optional, for symmetric proteins.

    /${protein-scaffold}_${ligand} // manually created
        /${protein-scaffold}-A_${ligand}-intermediate // Created by match.sh.
            ...
        /${protein-scaffold}-A_${ligand}-TS // Created by match.sh.
            ...
        /${protein-scaffold}-A_${ligand}-product // Created by match.sh.
            /${protein-scaffold}-A_${n} // n=1,2,3... Created by match.sh Contains Match output CloudPDB files. Will be removed after running generate_fast_design_input.sh.
            ...
            /match // Created by generate_fast_design_input.sh. Contains CloudPDB files for all matches moved from the ${protein-scaffold}-A_${n} folders.
            /X${position1}Z${position2} // Created by generate_fast_design_input.sh. Contains FastDesign input files of each match.
                X${position1}Z${position2}.pdb // Created by generate_fast_design_input.sh.
                ${ligand-3-letter-name}.rotlib.pdb // Created by generate_fast_design_input.sh.
                /design // Created by fast_design_matches.sh. Contains FastDesign output files for the match.
                ${protein-scaffold}-A_${point_mutation1}_${point_mutation2}_${point_mutation3}.pdb // Created by revert_designed_residues.sh when calling get_the_best_design.py. The best decoy of all FastDesign output decoys.
                /revert_${point_mutation3} // Revert each designed point mutations back to the wild type residue. Created by revert_designed_residues.sh.
            ...
        /${protein-scaffold}-AB_${ligand}-intermediate // created by match.sh
            ...
        /${protein-scaffold}-AB_${ligand}-TS // created by match.sh
            ...
        /${protein-scaffold}-AB_${ligand}-product // created by match.sh
            ...
            
