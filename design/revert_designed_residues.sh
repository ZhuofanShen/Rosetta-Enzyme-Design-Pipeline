#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";; # i.e. ../CPG2/CPG2-AB
        -linker) linker="$2";; # i.e. ../pAaF/pAaF-product
        -symm) symmetry="$2";; # optional
        -dup) duplicate_match="$2";; # optional
        -nbh) neighborhood="$2";; # optional
        -n) decoys="$2";; # optional
        -mem) memory="$2";;
        -ligand) ligand="$2";; # Optional, only for ligands or cofactors of the protein, not for matching residues or ligands.
        -not_get_best) not_get_the_best_design="$2";; # Optional. If this flag is given, get_the_best_design.py will not be called.
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
    decoys="-decoys 3"
else
    decoys="-decoys "${decoys}
fi

# if ! [ -z "${memory}" ]
# then
#     memory="--mem ${memory}"
# fi

path_to_protein=${position_files%/*}
scaffold=${position_files##*/}
substrate=${linker##*/}

if [ -z "${not_get_the_best_design}" ]
then
    linker_res_name=$(ls ${linker}/*_ligand.params)
    linker_res_name=${linker_res_name##*/}
    linker_res_name=${linker_res_name%_ligand.params}

    params_files="-params"
    for params_file in `ls ${linker}/*.params`
    do
        params_files=${params_files}" "${params_file}
    done

    IFS=','
    if ! [ -z "${ligand}" ]
    then
        params_files=${params_files}${ligand[@]}
    fi
    IFS='
    '

    protein=${scaffold%-*}
    reference_pdb=$(ls ${path_to_protein}/${protein}*.pdb)

    if ! [ -z "${symmetry}" ]
    then
        ignore_symmetric_mutations="-symm"
    fi

    # srun --job-name get_the_best_design --partition p_sdk94_1 --time 1:00:00 \
        python ../scripts/get_the_best_design.py ${scaffold}_${substrate} -ref ${reference_pdb} \
        -linker ${linker_res_name} ${params_files} ${ignore_symmetric_mutations}
fi

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

cd ${scaffold}_${substrate}

for variant in `ls`
do
    if [[ ${variant} == X*Z* ]]
    then
        cd ${variant}
        design=$(ls ${scaffold}_*.pdb)
        design=${design:0:-4}
        point_mutations_str=${design##${scaffold}_}
        IFS='_'
        point_mutations=(${point_mutations_str[@]})
        IFS='
        '
        for i in ${!point_mutations[@]}
        do
            if [[ ${point_mutations[${i}]} == *X ]] || [[ ${point_mutations[${i}]} == *Z ]]
            then
                unset 'point_mutations[${i}]'
            fi
        done
        if [[ ${#point_mutations[@]} > 0 ]]
        then
            for reverted_point_mutation in ${point_mutations[@]}
            do
                if [[ ${#point_mutations[@]} > 1 ]]
                then
                    point_mutations_argument="-muts"
                    for point_mutation in ${point_mutations[@]}
                    do
                        if [[ ${point_mutation} != ${reverted_point_mutation} ]]
                        then
                            point_mutations_argument=${point_mutations_argument}" "${point_mutation:1:-1}","${point_mutation: -1}
                        fi
                    done
                fi
                mkdir revert_${reverted_point_mutation}
                cd revert_${reverted_point_mutation}
                slurmit.py --job ${variant} --mem ${memory} --command \
                    "python ../../../../scripts/fast_design.py ../${variant}.pdb \
                    -sf ref2015_cst ${params_files} ${symmetry} \
                    -enzdescst ../../../${linker}/${substrate}${enzdes_cst_suffix} \
                    ${point_mutations_argument} ${neighborhood} -rmsd True ${decoys};"
                sleep 0.05
                cd ..
            done
        fi
        cd ..
    fi
done
cd ..
