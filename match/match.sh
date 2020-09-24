#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pdb) pdb="$2";;
        -ligand) ligand="$2";; #optional, only for ligands of the protein
        -positions) positions="$2";;
        -linker) linker="$2";;
        -mem) memory="$2";; #optional
        *) break;
    esac
    shift 2
done

if ! [ -z "${memory}" ]
then
  memory="--mem ${memory}"
fi

IFS=','
if ! [ -z "${ligand}" ]
then
  ligand="-extra_res_fa"
  for ligand_file in ${ligand[@]}
  do
      ligand=${ligand}" ../../"${ligand_file}
  done
fi
IFS='
'

scaffold=${positions##*/}
substrate=${linker##*/}
mkdir ${scaffold}_${substrate}
cd ${scaffold}_${substrate}

for pos in `ls ../${positions}/*.pos`
do
    scaffold_part=${pos##*/}
    scaffold_part=${scaffold_part%.pos}
    mkdir ${scaffold_part}
    cd ${scaffold_part}
    slurmit.py --job ${scaffold_part}_${substrate} ${memory} --command \
        "$ROSETTA3/bin/match.linuxgccrelease @../../../scripts/general_match.flags \
        @../../${linker}/subs.flags -match:scaffold_active_site_residues_for_geomcsts ../${pos} \
        -s ../../${pdb} ${ligand}"
    # pos file: ../${pos}
    # other input files: ../../${input}
    cd ..
    sleep 0.1
done