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
        mkdir match
        cp ${scaffold}*/* match
        # rm -rf ${scaffold}*
        # make sure that your python3 or anaconda3 has the pymol module installed.
        # conda config --set auto_activate_base false
        # eval "$(/home/zs251/anaconda3/bin/conda shell.bash hook)"
        python ../../scripts/generate_fast_design_input.py match ${homomeric}
        # conda deactivate
        cd ..
    fi
done

python ../scripts/find_match_intersection.py ${directories}
