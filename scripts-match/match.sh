#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";; # i.e. ../CPG2/CPG2-AB
        -lig) ligand="$2";; # i.e. ../pAaF/pAaF-product
        -mem) memory="$2";;
        *) break;
    esac
    shift 2
done

if ! [ -z "${memory}" ]
then
    memory="--mem ${memory}"
else
    memory="--mem 8000"
fi

path_to_protein=${position_files%/*}
IFS=','
params_files="-extra_res_fa "$(ls ${path_to_protein}/*.params)
if [[ ${params_files} == "-extra_res_fa " ]]
then
    params_files=""
fi
IFS='
'

scaffold=${position_files##*/}
protein=${scaffold%-*}
substrate=${ligand##*/}
reference_pdb=$(ls ${path_to_protein}/${protein}*.pdb)

mkdir ${scaffold}_${substrate}
cd ${scaffold}_${substrate}

for pos in `ls ../${position_files}/*.pos`
do
    scaffold_part=${pos##*/}
    scaffold_part=${scaffold_part%.pos}
    mkdir ${scaffold_part}
    cd ${scaffold_part}
    slurmit.py --job ${scaffold_part} ${memory} --command "~/Rosetta/main/source/bin/match.default.linuxgccrelease \
        @../../../../scripts-match/general_match.flags @../../${ligand}/subs.flags ${params_files} \
        -match:scaffold_active_site_residues_for_geomcsts ../${pos} -s ../../${reference_pdb}"
    cd ..
    sleep 0.05
done
