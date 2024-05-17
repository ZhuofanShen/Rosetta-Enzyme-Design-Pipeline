for pdb in `ls HEM`; do
    cd HEM/${pdb};
    for i in {1..1}; do
        mkdir relax$i;
        cd relax$i;
        slurmit.py --job ${pdb}${i} --command "python ../../../scripts/fast_design.py ../${pdb}_clean.pdb -ref ../${pdb}_clean.pdb -optH -sf ref2015_cst -no_min_sc_ids HEM -n 5 -prefix ${pdb} --score_terms fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1";
        cd ..;
    done
    cd ../..;
done
