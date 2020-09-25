#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";;
        -linker) linker="$2";; # i.e. ../pAaF/pAaF-product
        -symm) symmetry="$2";; #optional
        -nbh) neighborhood="$2";; #optional
        -n) decoys="$2";; #optional
        -mem) memory="$2";;
        -ligand) ligand="$2";; #optional, only for ligands or cofactors of the protein, not matching residues
        *) break;
    esac
    shift 2
done

if ! [ -z "${symmetry}" ]
then
    symmetry="-symm ../../../"${symmetry}
fi

if ! [ -z "${neighborhood}" ]
then
    neighborhood="-nbh "${neighborhood}
fi

if [ -z "${decoys}" ]
then
    decoys="-decoys 50"
else
    decoys="-decoys "${decoys}
fi

# if ! [ -z "${memory}" ]
# then
#     memory="--mem ${memory}"
# fi

linker_res_name=$(ls ${linker}/*_design.params)
linker_res_name=${linker_res_name##*/}
linker_res_name=${linker_res_name%_design.params}

IFS=','
if ! [ -z "${ligand}" ]
then
    params_files=""
    for params_file in ${ligand[@]}
    do
        params_files=${ligand}" ../../../"${params_file}
    done
fi
IFS='
'

scaffold=${position_files##*/}
substrate=${linker##*/}
cd ${scaffold}_${substrate}

for variant in `ls`
do
    if [[ ${variant} == X*Z* ]]
    then
        cd ${variant}
        mkdir design
        cd design
        slurmit.py --job ${variant} --mem ${memory} --command "python ../../../../scripts/fast_design.py \
            ../${variant}.pdb -sf beta_nov16_cst -params ../../../${linker}/${linker_res_name}_design.params \
            ../../../${linker}/CYX.params ../../../${linker}/TYZ.params ${params_files} ${symmetry} \
            -enzdescst ../../../${linker}/${substrate}_design.cst -rmsd True -nataa True \
            --enzdes ${neighborhood} ${decoys};"
        cd ../..
        sleep 0.1
    fi
done
cd ..
