#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";;
        -dup) duplicate_match="$2";; # optional
        -sym) symmetry="$2";; # optional
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

scaffold=${position_files##*/}
for scaffold_linker in `ls`
do
    if [[ ${scaffold_linker} == ${scaffold}_* ]]
    then
        directories=${directories}" "${scaffold_linker}
        cd ${scaffold_linker}
        for pos_part in `ls`
        do
            if [[ ${pos_part} == ${scaffold}_* ]]
            then
                python ../../scripts/generate_fast_design_input.py ${pos_part} ${duplicate_match} ${symmetry}
            fi
        done
        cd ..
    fi
done

python ../scripts/find_match_intersection.py ${directories}
