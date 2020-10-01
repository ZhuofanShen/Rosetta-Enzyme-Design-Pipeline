#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";; # i.e. ../CPG2/CPG2-AB
        -linker) linker="$2";; # i.e. ../pAaF/pAaF-product
        -mem) memory="$2";;
        -ligand) ligand="$2";; # Optional, only for ligands or cofactors of the protein, not for matching residues or ligands.
        *) break;
    esac
    shift 2
done

# if ! [ -z "${memory}" ]
# then
#   memory="--mem ${memory}"
# fi

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

path_to_protein=${position_files%/*}
scaffold=${position_files##*/}
protein=${scaffold%-*}
reference_pdb=$(ls ${path_to_protein}/${protein}*.pdb)

substrate=${linker##*/}
mkdir ${scaffold}_${substrate}
cd ${scaffold}_${substrate}

for pos in `ls ../${position_files}/*.pos`
do
    scaffold_part=${pos##*/}
    scaffold_part=${scaffold_part%.pos}
    mkdir ${scaffold_part}
    cd ${scaffold_part}
    slurmit.py --job ${scaffold_part} --mem ${memory} --command "~/Rosetta/main/source/bin/match.linuxgccrelease \
        @../../../scripts/general_match.flags @../../${linker}/subs.flags ${params_files} \
        -match:scaffold_active_site_residues_for_geomcsts ../${pos} -s ../../${reference_pdb}"
    cd ..
    sleep 0.1
done
