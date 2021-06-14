#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";;
        -dup) duplicate_match="$2";; # optional
        -symm) symmetry="$2";; # optional
        -debug) debug="$2";; # optional
        *) break;
    esac
    shift 2
done

if ! [ -z "${duplicate_match}" ]
then
    duplicate_match="--duplicate_match"
fi

if ! [ -z "${symmetry}" ]
then
    symmetry="--symmetry"
fi

if ! [ -z "${debug}" ]
then
    debug="--debug_mode"
fi

scaffold=${position_files##*/}
for scaffold_substrate in `ls`
do
    if [[ ${scaffold_substrate} == ${scaffold}_* ]]
    then
        substrate=${scaffold_substrate#*_}
        ligand=${substrate%-*}
        for ligand_params in `ls ../../ligands/${ligand}/${substrate}/*_ligand.params`
        do
            break
        done
        cd ${scaffold_substrate}
        for pos_part in `ls`
        do
            if [[ ${pos_part} == ${scaffold}_* ]]
            then
                python ../../../scripts-match/generate_fast_design_input.py ${pos_part} --ligand_params ../${ligand_params} ${duplicate_match} ${symmetry} ${debug}
            fi
        done
        cd ..
        directories=${directories}" "${scaffold_substrate}
    fi
done

python ../../scripts-match/find_match_intersection.py ${directories}
