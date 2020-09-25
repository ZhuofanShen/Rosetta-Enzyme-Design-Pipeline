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

scaffold=${positions##*/}
scaffold=${positions##*\\}
and="-and"
or="-or"
for scaffold_linker in `ls`
do
    if [[ ${scaffold_linker} == ${scaffold}* ]]
    then
        linker=${scaffold_linker##*_}
        if [[ ${linker} == *-*-* ]]
        then
            or=${or}" "${scaffold_linker}
        else
            and=${and}" "${scaffold_linker}
        fi
        cd ${scaffold_linker}
        mkdir match
        mv ${scaffold}*/* match
        rm -rf ${scaffold}*
        # make sure that your python3 or anaconda3 has the pymol module installed.
        # conda config --set auto_activate_base false
        # eval "$(~/anaconda3/bin/conda shell.bash hook)"
        python ../../scripts/generate_fast_design_input.py match ${homomeric}
        # conda deactivate
        cd ..
    fi
done

if [[ ${or} == "-or" ]]
then
    ../scripts/find_match_intersection.py ${and}
else
    ../scripts/find_match_intersection.py ${and} ${or}
fi
