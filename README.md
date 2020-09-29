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
        ${protein-scaffold}.symm // optional
        ${protein-scaffold}-A // position files for matching within chain A, created by generate_position_files.py
        ${protein-scaffold}-AB // position files for matching across chain A and chain B, created by generate_position_files.py
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
            ${ligand}-product_design.cst // optional, for symmetric proteins
