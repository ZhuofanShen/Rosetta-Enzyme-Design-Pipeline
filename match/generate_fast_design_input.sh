#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -pos) position_files="$2";;
        -homo) homomeric="$2";;
        *) break;
    esac
    shift 2
done

if ! [ -z "${homomeric}" ]
then
    homomeric="--homomeric"
fi

scaffold=${position_files##*/}
for scaffold_linker in `ls`
do
    if [[ ${scaffold_linker} == ${scaffold}* ]]
    then
        directories=${directories}" "${scaffold_linker}
        cd ${scaffold_linker}
        for pos_part in `ls`
        do
            if [[ ${pos_part} == ${scaffold}_* ]]
            then
                python ../../scripts/generate_fast_design_input.py ${pos_part} ${homomeric}
            fi
        done
        cd ..
    fi
done

python ../scripts/find_match_intersection.py ${directories}
