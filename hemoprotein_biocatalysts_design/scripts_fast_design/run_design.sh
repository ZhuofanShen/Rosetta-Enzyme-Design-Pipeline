for path in `ls HEM/*/cyclopropanation_styrene_EDA/*/*/*rot*/design.sh`; do cd ${path%/*}; bash ${path##*/}; sleep 0.01; cd ../../../../../..; done
