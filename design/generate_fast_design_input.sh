#!/bin/bash
while (( $# > 1 ))
do
    case $1 in
        -positions) positions="$2";;
        *) break;
    esac
    shift 2
done

scaffold=${positions##*/}
scaffold=${positions##*\\}
chain=${scaffold##*-}

if [[ ${#chain} > 1 ]]
then
    homomeric="--homomeric"
fi

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
        python3 ../../scripts/generate_fast_design_input.py match ${homomeric}
        cd ..
    fi
done

if [[ ${or} == "-or" ]]
then
    python3 ../scripts/find_match_intersection.py ${and}
else
    python3 ../scripts/find_match_intersection.py ${and} ${or}
fi