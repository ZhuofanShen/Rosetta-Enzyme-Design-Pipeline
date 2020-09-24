#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -scaf) scaffold="$2";; # i.e. CPG2-AB
        -linker) linker="$2";; # i.e. ../pAaF/pAaF-product
        -params) params="$2";; #optional, only for ligands of the protein
        -symm) symmetry="$2";; #optional
        -nbh) neighborhood="$2";; #optional
        -n) decoys="$2";; #optional
        -mem) memory="$2";; #optional
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

if ! [ -z "${memory}" ]
then
    memory="--mem ${memory}"
fi

IFS=','
if ! [ -z "${params}" ]
then
    ligand=""
    for params_file in ${params[@]}
    do
        ligand=${ligand}" ../../../"${params_file}
    done
fi
IFS='
'

substrate=${linker##*/}
cd ${scaffold}_${substrate}

for variant in `ls`
do
    if [[ ${variant} == X*Z* ]]
    then
        cd ${variant}
        rotlib=$(ls *.rotlib.pdb)
        linker_res_name=${rotlib%%.rotlib.pdb*}
        mkdir design
        cd design
        slurmit.py --job ${variant} ${memory} --command "python3 ../../../../scripts/fast_design_matches.py \
            ../${variant}.pdb -sf beta_nov16_cst -params ../../../${linker}/${linker_res_name}_design.params \
            ../../../${linker}/CYX.params ../../../${linker}/TYZ.params ${symmetry} ${duplicated_chains}\
            -enzdescst ../../../${linker}/${substrate}_design.cst -rmsd True -nataa True \
            --enzdes ${neighborhood} ${decoys};"
        cd ../..
        sleep 0.1
    fi
done
cd ..
