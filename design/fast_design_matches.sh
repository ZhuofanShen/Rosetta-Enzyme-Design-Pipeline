#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";; # i.e. ../CPG2/CPG2-AB
        -linker) linker="$2";; # i.e. ../pAaF/pAaF-product
        -cst_suffix) cst_suffix="$2";; # optional
        -nbh) neighborhood="$2";; #optional
        -n) decoys="$2";; #optional
        -mem) memory="$2";;
        *) break;
    esac
    shift 2
done

if ! [ -z "${cst_suffix}" ]
then
    enzdes_cst_suffix=_${cst_suffix}.cst
else
    enzdes_cst_suffix=.cst
fi

if ! [ -z "${neighborhood}" ]
then
    neighborhood="-nbh "${neighborhood}
fi

if [ -z "${decoys}" ]
then
    decoys="-n 50"
else
    decoys="-n "${decoys}
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
path_to_protein=${position_files%/*}
for params_file in `ls ${path_to_protein}/*.params`
do
    params_files=${params_files}" ../../../"${params_file}
done

for symmetry in `ls ${path_to_protein}/*.symm`
do
    symmetry="-symm ../../../"${symmetry}
    break
done

scaffold=${position_files##*/}
substrate=${linker##*/}
cd ${scaffold}_${substrate}

for variant in `ls`
do
    if [[ ${variant} == X*Z* ]] && [[ ${variant} != *_deprecated ]] && [ ! -f "${variant}/design" ]
    then
        cd ${variant}
        mkdir design
        cd design
        slurmit.py --job ${variant} --mem ${memory} --command \
            "python ../../../../scripts/fast_design.py ../${variant}.pdb \
            ${params_files} ${symmetry} -sf ref2015_cst \
            --score_terms fa_intra_rep_nonprotein:0.545 fa_intra_atr_nonprotein:1 \
            -enzdes_cst ../../../${linker}/${substrate}${enzdes_cst_suffix} \
            -enzdes -no_cys -nataa 1.5 ${neighborhood} ${decoys};"
        cd ../..
        sleep 0.05
    fi
done
cd ..
