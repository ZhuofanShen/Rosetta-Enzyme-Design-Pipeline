#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";; # i.e. ../CPG2/CPG2-AB
        -linker) linker="$2";; # i.e. ../pAaF/pAaF-product
        -symm) symmetry="$2";; #optional, symmetric file.
        -dup) duplicate_match="$2";; # optional
        -nbh) neighborhood="$2";; #optional
        -n) decoys="$2";; #optional
        -mem) memory="$2";;
        -ligand) ligand="$2";; # Optional, only for ligands or cofactors of the protein, not for matching residues or ligands.
        *) break;
    esac
    shift 2
done

if ! [ -z "${symmetry}" ]
then
    symmetry="-symm ../../../"${symmetry}
fi

if ! [ -z "${duplicate_match}" ]
then
    enzdes_cst_suffix=_symm.cst
else
    enzdes_cst_suffix=.cst
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

params_files="-params"
for params_file in `ls ${linker}/*.params`
do
    params_files=${params_files}" ../../../"${params_file}
done

IFS=','
if ! [ -z "${ligand}" ]
then
    for params_file in ${ligand[@]}
    do
        params_files=${params_files}" ../../../"${params_file}
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
            ../${variant}.pdb -sf ref2015_cst ${params_files} ${symmetry} \
            -enzdescst ../../../${linker}/${substrate}${enzdes_cst_suffix} \
            --enzdes -nataa True ${neighborhood} -rmsd True ${decoys};"
        cd ../..
        sleep 0.05
    fi
done
cd ..
