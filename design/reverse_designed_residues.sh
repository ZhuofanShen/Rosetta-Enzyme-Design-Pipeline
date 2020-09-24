#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pdb) reference_pdb="$2";;
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

pdb_name=
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
        design=$(ls ${scaffold}_*.pdb)
        design=${design:0:-4}
        design=${design##${scaffold}_}
        IFS='_'
        design=(${design[@]})
        IFS='
        '
        for i in ${!design[@]}
        do
            if [[ ${design[${i}]} == *X or ${design[${i}]} == *Z ]]
            then
                unset 'design[${i}]'
            fi
        done
        if [[ ${#design[@]} > 0 ]]
        then
            for reverted_point_mutation in ${design[@]}
            do
                if [[ ${#design[@]} > 1 ]]
                then
                    point_mutations="-muts"
                    for point_mutation in ${design[@]}
                    do
                        if [[ ${point_mutation} != ${reverted_point_mutation} ]]
                        then
                            point_mutations=${point_mutations}" "${point_mutation:1:-1}","${point_mutation: -1}
                        fi
                    done
                else
                    point_mutations=""
                fi
                mkdir revert_${reverted_point_mutation}
                cd revert_${reverted_point_mutation}
                slurmit.py --job ${variant} ${memory} --command "python3 ../../../../scripts/fast_design_matches.py \
                    ../${variant}.pdb -sf beta_nov16_cst -params ../../../${linker}/${linker_res_name}_design.params \
                    ../../../${linker}/CYX.params ../../../${linker}/TYZ.params ${symmetry} ${duplicated_chains}\
                    -enzdescst ../../../${linker}/${substrate}_design.cst -rmsd True -nataa True \
                    ${point_mutations} ${neighborhood} ${decoys};"
                sleep 0.1
                cd ..
            done
        fi
        cd ..
    fi
done
cd ..