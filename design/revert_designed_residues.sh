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
        -get_best) get_the_best_design="$2";; # Optional. If this flag is given, get_the_best_design.py will be called.
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
        params_files=${ligand}" ../"${params_file}
    done
fi

path_to_protein=${position_files%/*}
scaffold=${position_files##*/}
substrate=${linker##*/}

if ! [ -z "${get_the_best_design}" ]
then
    protein=${scaffold%-*}
    reference_pdb=$(ls ${path_to_protein}/${protein}*.pdb)

    if ! [ -z "${symmetry}" ]
    then
        homomeric="-homo"
    fi

    srun --job-name get_the_best_design --partition p_sdk94_1 --time 1:00:00 \
        python ../scripts/get_the_best_design.py ${scaffold}_${substrate} -ref ${reference_pdb} \
        -params ${linker}/${linker_res_name}_design.params ${linker}/CYX.params \
        ${linker}/TYZ.params ${params_files} -linker ${linker_res_name} ${homomeric}
fi

if ! [ -z "${ligand}" ]
then
    params_files=""
    for params_file in ${ligand[@]}
    do
        params_files=${ligand}" ../../"${params_file}
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
        design=${design##${scaffold}_}
        IFS='_'
        design=(${design[@]})
        IFS='
        '
        for i in ${!design[@]}
        do
            if [[ ${design[${i}]} == *X ]] || [[ ${design[${i}]} == *Z ]]
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
                slurmit.py --job ${variant} --mem ${memory} --command "python ../../../../scripts/fast_design.py \
                    ../${variant}.pdb -sf ref2015_cst -params ../../../${linker}/${linker_res_name}_design.params \
                    ../../../${linker}/CYX.params ../../../${linker}/TYZ.params ${params_files} ${symmetry} \
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
